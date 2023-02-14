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
  Assert(_distinct_count_per_bin.size() == _bin_heights.size(), "Must have the same number of distinct counts and heights.");

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
      -------------------------------------   
      |  Sets current minima and maxima.  |
      -------------------------------------
    */
    std::vector<T> bin_minima(barrier_count + 1);
    std::vector<T> bin_maxima(barrier_count + 1);
    bin_minima[0] = value_distribution[0].first;
    bin_maxima[bin_maxima.size() - 1] = value_distribution[value_distribution.size() - 1].first;
    { // barrier_ind has to count outside each iteration, but should ONLY be used inside the loop. 
      // Therefore, we use a block.
      auto barrier_ind = uint32_t{0};
      for (auto val_ind = uint32_t{0}; val_ind < value_distribution.size() - 2; ++val_ind) {
          if (val_ind == barrier_indexes[barrier_ind]) {
              bin_maxima[barrier_ind] = value_distribution[val_ind].first;
              bin_minima[barrier_ind + 1] = value_distribution[val_ind + 1].first;
              ++barrier_ind;
          }
      }
    }

    /*
      -----------------------
      |  Error Calculation  |
      -----------------------
    */
    for (auto barrier_index = uint32_t{0}; barrier_index < barrier_count; ++barrier_index) {

        /*
          -----------------------------------------------------
          |  Calculate error increase for removing a barrier  |
          -----------------------------------------------------
        */
        error_increase ei;
        ei.barrier_index = barrier_index;
        {
          const auto current_error = static_cast<float>(bin_maxima[barrier_index] - bin_minima[barrier_index]);
          const auto current_error_of_next_bin = static_cast<float>(
              bin_maxima[barrier_index + 1] - bin_minima[barrier_index + 1]
          );

          const auto total_current_error = std::pow(std::max(current_error, current_error_of_next_bin), 2);
          auto new_error = std::pow(static_cast<float>(bin_maxima[barrier_index + 1] - bin_minima[barrier_index]), 2);

          // Error increase: Error of new bin (if barrier removed) - error of old bin (if barrier not removed).
          ei.error_increase = static_cast<float>(std::sqrt(new_error - total_current_error));
          error_increases.emplace_back(ei);
        }

        /*
          -----------------------------------------------------
          |  Calculate error decrease for removing a barrier  |
          -----------------------------------------------------
        */
        error_decrease ed;
        ed.ideal_barrier_index = value_distribution.size() + 1; // Set to invalid value.
        ed.error_decrease = 0;
        {
          const auto barrier_position = barrier_indexes[barrier_index];
          // Find next barrier. Might also be the last value, if there is no next barrier.
          auto next_barrier_position = ((barrier_position + 1) < barrier_indexes.size()) 
                                        ? barrier_indexes[barrier_index + 1] 
                                        : value_distribution.size() - 1;

          // Start looking at the value after the current barrier.
          auto val_position = barrier_position + 1;

          while (val_position < next_barrier_position) {
              // Calculate for each index, how the error would be influenced.
              const auto& [value, frequency] = value_distribution[val_position];

              const auto curr_error = static_cast<float>(bin_maxima[barrier_index+1] - bin_minima[barrier_index+1]);
              const auto new_left_error = static_cast<float>(value - bin_minima[barrier_index+1]);
              const auto new_right_error = static_cast<float>(bin_maxima[barrier_index + 1] - value);

              const auto error_decrease = static_cast<float>(std::sqrt(
                  std::pow(curr_error, 2) - std::pow(std::max(new_left_error, new_right_error), 2)
              ));

              if (error_decrease > ed.error_decrease) {
                  ed.error_decrease = error_decrease;
                  ed.ideal_barrier_index = val_position;
              }
              ++val_position;
          }
        }

        if (ed.ideal_barrier_index != value_distribution.size() + 1 && ed.error_decrease > 0) {
            error_decreases.emplace_back(ed);
        }

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

  auto binCount = max_bin_count;
  if (max_bin_count > value_distribution.size()) {
    binCount = value_distribution.size();
  }

  // Trivial Histogram.
  if (value_distribution.size() == 1) {
      std::vector<T> bin_minima(1);
      bin_minima[0] = value_distribution[0].first;
      std::vector<T> bin_maxima(1);
      bin_maxima[0] = value_distribution[0].first;
      std::vector<HistogramCountType> bin_heights(1);
      bin_heights[0] = value_distribution[0].second;
      std::vector<HistogramCountType> bin_distinct_counts(1);
      bin_distinct_counts[0] = 1;

      return std::make_shared<GDYHistogram<T>>(
              std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
              std::move(bin_distinct_counts), total_count, domain);
  }


  std::vector<uint32_t> barrier_indexes(binCount - 1);
  // Removing barriers.
  std::vector<error_increase> error_increases(binCount - 1);
  // Adding barriers in segments.
  std::vector<error_decrease> error_decreases(binCount - 1);

  // 1. Generate random histogram (aka value distribution partitioning).
  //auto histogram = EqualDistinctCountHistogram<T>::from_column(table, column_id, max_bin_count, domain);
  for (auto i = 0u; i < binCount - 1; ++i) {
    barrier_indexes[i] = i;
  }

  // 2. Calculate errors for removing barriers and for adding an ideal partitioning.
  auto populated_errors = calculate_error_changes(value_distribution, barrier_indexes);

  // These are sorted. Means, the largest error increase & decrease is at the beginning.
  error_increases = populated_errors.first;
  error_decreases = populated_errors.second;

  // --------------------------------------------------------
  // -------          Construct Histogram             -------
  // --------------------------------------------------------

  auto retries = uint32_t{0};
  while (error_decreases.size() > 0 && error_increases[0].error_increase < error_decreases[0].error_decrease && static_cast<float>(retries) < total_count) {
    auto largest_error_increase = error_increases[0];
    auto largest_error_decrease = error_decreases[0];

    Assert(barrier_indexes[largest_error_increase.barrier_index] != largest_error_decrease.ideal_barrier_index, "Barrier cannot be repositioned to itself.");

    // Replace barrier with better one.
    barrier_indexes[largest_error_increase.barrier_index] = largest_error_decrease.ideal_barrier_index;
    // Sort barrier indexes, since after the last exchange, they are not anymore in order.
    std::sort(barrier_indexes.begin(), barrier_indexes.end());

    populated_errors = calculate_error_changes(value_distribution, barrier_indexes);
    error_increases = populated_errors.first;
    error_decreases = populated_errors.second;

    ++retries;
  }

  std::vector<T> bin_minima(binCount);
  std::vector<T> bin_maxima(binCount);
  std::vector<HistogramCountType> bin_heights(binCount);
  std::vector<HistogramCountType> bin_distinct_counts(binCount);

  // Populate bin_minima & bin_maxima.
  auto bin_id = BinID{0};
  bin_minima[bin_id] = value_distribution[0].first;
  for (auto val_ind = uint32_t{0}; val_ind < value_distribution.size() - 1; ++val_ind) {
    if (val_ind == barrier_indexes[bin_id]) {
      bin_maxima[bin_id] = value_distribution[val_ind].first;
      ++bin_id;
      bin_minima[bin_id] = value_distribution[val_ind + 1].first;
    }
    bin_heights[bin_id] += value_distribution[val_ind].second;
    ++bin_distinct_counts[bin_id];
  }
  // Last value is not yet counted in.
  bin_maxima[binCount - 1] = value_distribution.back().first;
  bin_heights[binCount - 1] += value_distribution.back().second;
  ++bin_distinct_counts[binCount - 1];

  return std::make_shared<GDYHistogram<T>>(
      std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
      std::move(bin_distinct_counts), total_count, domain);
}

template <typename T>
std::string GDYHistogram<T>::name() const {
  return "GDYHistogram";
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> GDYHistogram<T>::clone() const {
    // The new histogram needs a copy of the data
    auto bin_minima_copy = _bin_minima;
    auto bin_maxima_copy = _bin_maxima;
    auto bin_heights_copy = _bin_heights;
    auto bin_distinct_counts_copy = _distinct_count_per_bin;

    Assert(_bin_minima.size() == _bin_maxima.size(), "CLONING: Must have the same number of lower as upper bin edges.");
    Assert(_bin_minima.size() == _bin_heights.size(), "CLONING: Must have the same number of edges and heights.");
    Assert(_distinct_count_per_bin.size() == _bin_heights.size(), "CLONING: Must have the same number of distinct counts and heights.");

    return std::make_shared<GDYHistogram<T>>(std::move(bin_minima_copy), std::move(bin_maxima_copy),
                                             std::move(bin_heights_copy), std::move(bin_distinct_counts_copy),
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
  if (bin_count() == 1) {
    return 1;
  }

  return _distinct_count_per_bin.at(index);
}

template <typename T>
HistogramCountType GDYHistogram<T>::total_count() const {
  return _total_count;
}

template <typename T>
HistogramCountType GDYHistogram<T>::total_distinct_count() const {
  return _total_distinct_count;
}

template class GDYHistogram<int32_t>;
template class GDYHistogram<int64_t>;
template class GDYHistogram<float>;
template class GDYHistogram<double>;

}  // namespace hyrise
