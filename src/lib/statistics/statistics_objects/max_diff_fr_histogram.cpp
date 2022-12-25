#include "max_diff_fr_histogram.hpp"

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

namespace {

using namespace hyrise;  // NOLINT

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
MaxDiffFrHistogram<T>::MaxDiffFrHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                                            std::vector<HistogramCountType>&& bin_height, std::vector<HistogramCountType>&& bin_distinct_counts,
                                            const HistogramCountType total_count, const HistogramDomain<T>& domain)
    : AbstractHistogram<T>(domain),
    _bin_minima(std::move(bin_minima)),
      _bin_maxima(std::move(bin_maxima)),
      _bin_heights(std::move(bin_height)),
      _bin_distinct_counts(std::move(bin_distinct_counts)),
      _total_count{total_count} {}

template <typename T>
std::string MaxDiffFrHistogram<T>::name() const {
  return "MaxDiffFr";
}

template <typename T>
std::shared_ptr<MaxDiffFrHistogram<T>> MaxDiffFrHistogram<T>::from_column(const Table& table,
                                                                            const ColumnID column_id,
                                                                            const BinID max_bin_count,
                                                                            const HistogramDomain<T>& domain) {
  Assert(max_bin_count > 0, "max_bin_count must be greater than zero ");

  // Compute the value distribution. Basically, counting how many times each value appears in
  // the column.
  const auto value_distribution = value_distribution_from_column(table, column_id, domain);

  if (value_distribution.empty()) {
    return nullptr;
  }

  // Get the total number of values present in the Histogram.
  auto total_count = HistogramCountType{0};
  for (const auto& value_freq : value_distribution) {
    total_count += value_freq.second;
  }

  // Find the number of bins.
  auto bin_count = max_bin_count;
  if (static_cast<BinID>(value_distribution.size()) < max_bin_count) {
    bin_count = static_cast<BinID>(value_distribution.size());
  }

  // Initialize the resulting data structures.
  std::vector<T> bin_minima(bin_count);
  std::vector<T> bin_maxima(bin_count);
  std::vector<HistogramCountType> bin_heights(bin_count, 0);
  std::vector<HistogramCountType> bin_distinct_counts(bin_count, 0);

  // TODO(robert, tobias): What is the best beta?
  // Currently, the average frequency of a term is used.
  const auto beta = std::ceil(total_count / static_cast<float>(value_distribution.size()));

  auto bin_index = BinID{0};
  bin_minima[0] = value_distribution[0].first;
  bin_maxima[0] = value_distribution[0].first;
  bin_heights[0] = value_distribution[0].second;
  bin_distinct_counts[0] = 1;

  for (auto ind = uint32_t{1}; ind < value_distribution.size(); ++ind) {
    const auto val = value_distribution[ind];
    auto prevVal = value_distribution[ind - 1];

    // If the "difference" between two values is larger than beta --> new bin.
    const auto distance = std::abs(val.second - prevVal.second);
    if (distance >= beta && bin_index < bin_count) {
        // Next Bin
        ++bin_index;
        bin_minima[bin_index] = val.first;
    }

    bin_maxima[bin_index] = val.first;
    bin_heights[bin_index] += val.second;
    bin_distinct_counts[bin_index] += 1;
  }

  auto usedBuckets = 0;
  for (const auto height : bin_heights) {
    if (height > 0) ++usedBuckets;
  }

  // Initialize the resulting data structures.
  std::vector<T> actual_bin_minima(usedBuckets);
  std::vector<T> actual_bin_maxima(usedBuckets);
  std::vector<HistogramCountType> actual_bin_heights(usedBuckets, 0);
  std::vector<HistogramCountType> actual_bin_distinct_counts(usedBuckets, 0);
  for (auto ind = 0; ind < usedBuckets; ++ind) {
    actual_bin_minima[ind] = bin_minima[ind];
    actual_bin_maxima[ind] = bin_maxima[ind];
    actual_bin_heights[ind] = bin_heights[ind];
    actual_bin_distinct_counts[ind] = bin_distinct_counts[ind];
  }

  return std::make_shared<MaxDiffFrHistogram<T>>(std::move(actual_bin_minima), std::move(actual_bin_maxima), std::move(actual_bin_heights), std::move(actual_bin_distinct_counts), total_count);
}

template <typename T>
HistogramCountType MaxDiffFrHistogram<T>::total_count() const {
  return _total_count;
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> MaxDiffFrHistogram<T>::clone() const {
  // The new histogram needs a copy of the data
  auto bin_minima_copy = _bin_minima;
  auto bin_maxima_copy = _bin_maxima;
  auto bin_heights_copy = _bin_heights;
  auto bin_distinct_counts_copy = _bin_distinct_counts;

  return std::make_shared<MaxDiffFrHistogram<T>>(std::move(bin_minima_copy), std::move(bin_maxima_copy),
                                                  std::move(bin_heights_copy), std::move(bin_distinct_counts_copy), _total_count);
}

template <typename T>
BinID MaxDiffFrHistogram<T>::bin_count() const {
  return _bin_maxima.size();
}

template <typename T>
const T& MaxDiffFrHistogram<T>::bin_minimum(BinID index) const {
  return _bin_minima[index];
}

template <typename T>
const T& MaxDiffFrHistogram<T>::bin_maximum(BinID index) const {
  return _bin_maxima[index];
}

template <typename T>
HistogramCountType MaxDiffFrHistogram<T>::bin_height(const BinID index) const {
  return _bin_heights[index];
}

template <typename T>
BinID MaxDiffFrHistogram<T>::_bin_for_value(const T& value) const {
  for (auto index = size_t{0}; index < this->bin_count() - 1; ++index) {
    if (bin_minimum(index + 1) > value) {
      return index;
    }
  }
  // Value must be in last bin.
  return this->bin_count() - 1;
}

template <typename T>
BinID MaxDiffFrHistogram<T>::_next_bin_for_value(const T& value) const {
  const auto bin_for_value = this->_bin_for_value(value);
  // If value is in last bin, return first bin.
  if (bin_for_value == bin_count() - 1) {
    return 0;
  }
  return bin_for_value + 1;
}

template <typename T>
HistogramCountType MaxDiffFrHistogram<T>::total_distinct_count() const {
  return _total_count;
}

template <typename T>
HistogramCountType MaxDiffFrHistogram<T>::bin_distinct_count(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");
  return _bin_distinct_counts[index];
}

EXPLICITLY_INSTANTIATE_DATA_TYPES(MaxDiffFrHistogram);

}  // namespace hyrise
