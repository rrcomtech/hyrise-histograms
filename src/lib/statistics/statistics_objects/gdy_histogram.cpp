#include "gdy_histogram.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "boost/sort/sort.hpp"
#include "tsl/robin_map.h"

#include "generic_histogram.hpp"
#include "resolve_type.hpp"
#include "storage/segment_iterate.hpp"
#include "statistics/statistics_objects/equal_distinct_count_histogram.hpp"

namespace {

using namespace hyrise;  // NOLINT

// Think of this as an unordered_map<T, HistogramCountType>. The hash, equals, and allocator template parameter are
// defaults so that we can set the last parameter. It controls whether the hash for a value should be cached. Doing
// so reduces the cost of rehashing at the cost of slightly higher memory consumption. We only do it for strings,
// where hashing is somewhat expensive.
template <typename T>
using ValueDistributionMap =
    tsl::robin_map<T, HistogramCountType, std::hash<T>, std::equal_to<T>,
                   std::allocator<std::pair<T, HistogramCountType>>, std::is_same_v<std::decay_t<T>, pmr_string>>;

template <typename T>
void add_segment_to_value_distribution(const AbstractSegment& segment, ValueDistributionMap<T>& value_distribution,
                                       const HistogramDomain<T>& domain) {
  segment_iterate<T>(segment, [&](const auto& iterator_value) {
    if (iterator_value.is_null()) {
      return;
    }

    if constexpr (std::is_same_v<T, pmr_string>) {
      // Do "contains()" check first to avoid the string copy incurred by string_to_domain() where possible
      if (domain.contains(iterator_value.value())) {
        ++value_distribution[iterator_value.value()];
      } else {
        ++value_distribution[domain.string_to_domain(iterator_value.value())];
      }
    } else {
      ++value_distribution[iterator_value.value()];
    }
  });
}

template <typename T>
std::vector<std::pair<T, HistogramCountType>> value_distribution_from_column(const Table& table,
                                                                             const ColumnID column_id,
                                                                             const HistogramDomain<T>& domain) {
  auto value_distribution_map = ValueDistributionMap<T>{};
  const auto chunk_count = table.chunk_count();
  for (auto chunk_id = ChunkID{0}; chunk_id < chunk_count; ++chunk_id) {
    const auto chunk = table.get_chunk(chunk_id);
    if (!chunk) {
      continue;
    }

    add_segment_to_value_distribution<T>(*chunk->get_segment(column_id), value_distribution_map, domain);
  }

  auto value_distribution =
      std::vector<std::pair<T, HistogramCountType>>{value_distribution_map.begin(), value_distribution_map.end()};
  value_distribution_map.clear();  // Maps can be large and sorting slow. Free space early.
  boost::sort::pdqsort(value_distribution.begin(), value_distribution.end(),
                       [&](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

  return value_distribution;
}
}  // namespace

namespace hyrise {

template <typename T>
GDYHistogram<T>::GDYHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                                                            std::vector<HistogramCountType>&& bin_heights,
                                                            std::vector<HistogramCountType>&&distinct_count_per_bin,
                                                            const BinID bin_count_with_extra_value,
                                                            const HistogramDomain<T>& domain)
    : AbstractHistogram<T>(domain),
      _bin_minima(std::move(bin_minima)),
      _bin_maxima(std::move(bin_maxima)),
      _bin_heights(std::move(bin_heights)),
      _distinct_count_per_bin(std::move(distinct_count_per_bin)),
      _bin_count_with_extra_value(bin_count_with_extra_value),
      _domain(domain) {
  Assert(_bin_minima.size() == _bin_maxima.size(), "Must have the same number of lower as upper bin edges.");
  Assert(_bin_minima.size() == _bin_heights.size(), "Must have the same number of edges and heights.");

  AbstractHistogram<T>::_assert_bin_validity();

  _total_count = std::accumulate(_bin_heights.cbegin(), _bin_heights.cend(), HistogramCountType{0});
}

/**
 * Calculates for each barrier the error according to GDY partitioning.
 *
 * Idea:
 *      Move a barrier to a position, so that the error is decreased as much as possible.
 *      Therefore, we need to know which barrier removal would lead to the smallest error increase
 *      and which barrier insertion would lead to the largest error decrease.
 *      Once calculated, we can move the barrier accordingly. This is done as long as there would
 *      be an error decrease.
 *
 * Note:
 *      This algorithm is a greedy algorithm. It does not guarantee to find the optimal solution.
 *      Also, an initial partitioning is required. This is done by the caller. Here, one can use
 *      other histograms (e.g. EquiHeight, EquiWidth, EqualDistinctCount) or a random partitioning.
 *      Author of GDY Paper: "The greedy boundary modiﬁcation operations matter for the ﬁnal result
 *      more than the choice of an intial partitioning" (Halim et al., 01/2009).
 *
 * @tparam T
 *          The type of the histogram domain. TODO: Also support strings.
 * @param value_distribution
 *          The value distribution of the column. Pairs of type (value, count).
 *          Access via value_distribution[i].first and value_distribution[i].second.
 * @param barrier_indexes
 *          The predefined barrier indexed. A barrier at index i means a barrier
 *          between value_distribution[i] and value_distribution[i+1].
 * @return
 *          A pair of vectors. The first vector contains the error increase for each barrier
 *          in case they would be removed.
 *          The second vector contains the error decrease and the ideal barriers where to insert
 *          a new barrier.
 */
template <typename T>
std::pair<std::vector<error_increase>, std::vector<error_decrease>> calculate_error_changes(
        std::vector<std::pair<T, HistogramCountType>> value_distribution, std::vector<uint32_t> barrier_indexes
    ) {
    const auto barrier_count = barrier_indexes.size();

    std::vector<error_increase> error_increases(0);
    std::vector<error_decrease> error_decreases(0);

    /*
     *  Set current minima and maxima.
     */
    std::vector<T> bin_minima(barrier_count + 1);
    std::vector<T> bin_maxima(barrier_count + 1);
    bin_minima[0] = value_distribution[0].first;
    auto barrier_ind = uint32_t{0};
    for (auto val_ind = uint32_t{0}; val_ind < value_distribution.size() - 1; ++val_ind) {
        bin_maxima[barrier_ind] = value_distribution[val_ind].first;
        if (barrier_indexes[barrier_ind] == val_ind) {
            ++barrier_ind;
            bin_minima[barrier_ind] = value_distribution[val_ind].first;
        }
    }

    for (auto ind = uint32_t{0}; ind < barrier_count; ++ind) {
        const auto curr_barrier = barrier_indexes[barrier_ind];
        const auto current_error = static_cast<float>(std::abs(bin_maxima[ind] - bin_minima[ind]));
        const auto current_error_of_next_bin = static_cast<float>(
            std::abs(bin_maxima[ind + 1] - bin_minima[ind + 1])
        );
        /*
         * The error increase can be viewed in comparison to two different bins.
         * It feels more natural to use the smaller error increase.
         */
        const auto total_current_error = std::pow(std::max(current_error, current_error_of_next_bin), 2);
        auto new_error = std::pow(static_cast<float>(std::abs(bin_maxima[ind + 1] - bin_minima[ind])), 2);
        error_increase ei;
        ei.barrier_index = curr_barrier;
        // Error increase: Error of new bin (if barrier removed) - error of old bin (if barrier not removed).
        ei.error_increase = static_cast<float>(new_error - total_current_error);

        // Evaluate ideal positioning for new barrier.
        auto val_partition_index = curr_barrier;

        auto next_barrier_index = uint32_t{0};
        if (curr_barrier + 1 < barrier_indexes.size()) {
            next_barrier_index = barrier_indexes[curr_barrier + 1];
        } else {
            next_barrier_index = value_distribution.size() - 1;
        }

        error_decrease ed;
        ed.ideal_barrier_index = 0;
        ed.error_decrease = 0;
        while ((val_partition_index+1) < next_barrier_index) {
            // Calculate for each index, how the error would be influenced.
            const auto& [partition_value, partition_count] = value_distribution[val_partition_index];

            const auto curr_error = static_cast<float>(std::abs(bin_maxima[ind + 1] - bin_minima[ind]));
            const auto new_left_error = static_cast<float>(std::abs(partition_value - bin_minima[ind]));
            const auto new_right_error = static_cast<float>(std::abs(bin_maxima[ind + 1] - partition_value));

            const auto error_decrease = static_cast<float>(
                std::pow(curr_error, 2) - (std::pow(new_left_error, 2) + std::pow(new_right_error, 2))
            );

            if (error_decrease > ed.error_decrease) {
                ed.error_decrease = error_decrease;
                ed.ideal_barrier_index = val_partition_index;
            }

            ++val_partition_index;
        }

        auto ei_contained = false;
        auto ed_contained = false;
        for (auto error_inc : error_increases) {
            if (error_inc.barrier_index == ei.barrier_index) {
                ei_contained = true;
                break;
            }
        }
        for (auto error_dec : error_decreases) {
            if (error_dec.ideal_barrier_index == ed.ideal_barrier_index) {
                ed_contained = true;
                break;
            }
        }
        if (!ei_contained) error_increases.emplace_back(ei);
        if (!ed_contained) error_decreases.emplace_back(ed);
    }

    // Sort error increases and decreases by error increase.
    std::sort(error_increases.begin(), error_increases.end(),
          [](error_increase ei1, error_increase ei2) {
              return ei1.error_increase < ei2.error_increase;
          });
    std::sort(error_decreases.begin(), error_decreases.end(),
          [](error_decrease ed1, error_decrease ed2) {
              return ed1.error_decrease > ed2.error_decrease;
          });

    return {error_increases, error_decreases};
}

template <typename T>
std::shared_ptr<GDYHistogram<T>> GDYHistogram<T>::from_column(
    const Table& table, const ColumnID column_id, const BinID max_bin_count, const HistogramDomain<T>& domain) {
  Assert(max_bin_count > 0, "max_bin_count must be greater than zero ");

  auto value_distribution = value_distribution_from_column(table, column_id, domain);
  auto total_count = float{0};
  for (const auto& [value, count] : value_distribution) {
    total_count += count;
  }

  if (value_distribution.empty()) {
    return nullptr;
  }

  std::vector<uint32_t> barrier_indexes(max_bin_count - 1);
  // Removing barriers.
  std::vector<error_increase> error_increases(max_bin_count - 1);
  // Adding barriers in segments.
  std::vector<error_decrease> error_decreases(max_bin_count - 1);

  // 1. Generate random histogram (aka value distribution partitioning).
  //auto histogram = EqualDistinctCountHistogram<T>::from_column(table, column_id, max_bin_count, domain);

  // As an initial partitioning, we use the values 1 to max_bin_count.
  const auto initialHistogram = EqualDistinctCountHistogram<T>::from_column(table, column_id, max_bin_count, domain);

  for (auto i = 0u; i < max_bin_count - 1; ++i) {
    barrier_indexes[i] = initialHistogram->bin_minimum(i + 1);
  }

  // 2. Calculate errors for removing barriers and for adding an ideal partitioning.
  auto populated_errors = calculate_error_changes(value_distribution, barrier_indexes);

  // These are sorted. Means, the largest error increase & decrease is at the beginning.
  error_increases = populated_errors.first;
  error_decreases = populated_errors.second;

  // --------------------------------------------------------
  // -------          Construct Histogram             -------
  // --------------------------------------------------------

  /*
   * This was written by GH Copilot. The student, who was supposed to write this, was too lazy
   * (the last sentence was suggested by Copilot ...).
   */
  while (error_increases[0].error_increase < error_decreases[0].error_decrease) {
    // 3.1. Remove barrier with largest error increase.
    auto largest_error_increase = error_increases[0];
    barrier_indexes.erase(barrier_indexes.begin() + largest_error_increase.barrier_index);
    // 3.2. Add barrier with largest error decrease.
    auto largest_error_decrease = error_decreases[0];
    auto largest_error_decrease_index = largest_error_decrease.ideal_barrier_index;
    barrier_indexes.insert(barrier_indexes.begin() + largest_error_decrease_index, largest_error_decrease_index);
    // 3.3. Recalculate errors.
    populated_errors = calculate_error_changes(value_distribution, barrier_indexes);
    error_increases = populated_errors.first;
    error_decreases = populated_errors.second;
  }

  std::vector<T> bin_minima(max_bin_count);
  std::vector<T> bin_maxima(max_bin_count);
  std::vector<HistogramCountType> bin_heights(max_bin_count);
  std::vector<HistogramCountType> bin_distinct_counts(max_bin_count);

  // Populate bin_minima & bin_maxima.
  auto bin_id = BinID{0};
  for (auto val_ind = uint32_t{0}; val_ind < value_distribution.size(); ++val_ind) {
    if (val_ind == barrier_indexes[bin_id]) {
      bin_maxima[bin_id] = value_distribution[val_ind].first;
      bin_id++;
      bin_minima[bin_id] = value_distribution[val_ind].first;
    }
    bin_heights[bin_id] += value_distribution[val_ind].second;
    ++bin_distinct_counts[bin_id];
  }
  bin_maxima[max_bin_count - 1] = value_distribution.back().first;

  return std::make_shared<GDYHistogram<T>>(
      std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
      std::move(bin_distinct_counts), total_count, domain);
}

template <typename T>
std::string GDYHistogram<T>::name() const {
  return "EqualDistinctCount";
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> GDYHistogram<T>::clone() const {
    // The new histogram needs a copy of the data
    auto bin_minima_copy = _bin_minima;
    auto bin_maxima_copy = _bin_maxima;
    auto bin_heights_copy = _bin_heights;

    return std::make_shared<GDYHistogram<T>>(std::move(bin_minima_copy), std::move(bin_maxima_copy),
                                             std::move(bin_heights_copy), std::move(_distinct_count_per_bin),
                                             total_count(), _domain);
}

template <typename T>
BinID GDYHistogram<T>::bin_count() const {
  return _bin_heights.size();
}

template <typename T>
BinID GDYHistogram<T>::_bin_for_value(const T& value) const {
  const auto iter = std::lower_bound(_bin_maxima.cbegin(), _bin_maxima.cend(), value);
  const auto index = static_cast<BinID>(std::distance(_bin_maxima.cbegin(), iter));

  if (iter == _bin_maxima.cend() || value < bin_minimum(index) || value > bin_maximum(index)) {
    return INVALID_BIN_ID;
  }

  return index;
}

template <typename T>
BinID GDYHistogram<T>::_next_bin_for_value(const T& value) const {
  const auto it = std::upper_bound(_bin_maxima.cbegin(), _bin_maxima.cend(), value);

  if (it == _bin_maxima.cend()) {
    return INVALID_BIN_ID;
  }

  return static_cast<BinID>(std::distance(_bin_maxima.cbegin(), it));
}

template <typename T>
const T& GDYHistogram<T>::bin_minimum(const BinID index) const {
  DebugAssert(index < _bin_minima.size(), "Index is not a valid bin.");
  return _bin_minima[index];
}

template <typename T>
const T& GDYHistogram<T>::bin_maximum(const BinID index) const {
  DebugAssert(index < _bin_maxima.size(), "Index is not a valid bin.");
  return _bin_maxima[index];
}

template <typename T>
HistogramCountType GDYHistogram<T>::bin_height(const BinID index) const {
  DebugAssert(index < _bin_heights.size(), "Index is not a valid bin.");
  return _bin_heights[index];
}

template <typename T>
HistogramCountType GDYHistogram<T>::bin_distinct_count(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");
  return _distinct_count_per_bin[index];
}

template <typename T>
HistogramCountType GDYHistogram<T>::total_count() const {
  return _total_count;
}

template <typename T>
HistogramCountType GDYHistogram<T>::total_distinct_count() const {
  return _total_distinct_count;
}

template class GDYHistogram<int>;
template class GDYHistogram<long>;
template class GDYHistogram<float>;
template class GDYHistogram<double>;

}  // namespace hyrise
