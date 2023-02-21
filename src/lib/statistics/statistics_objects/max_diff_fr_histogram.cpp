#include "max_diff_fr_histogram.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <fstream>

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
      _total_count{total_count} {

        for (const auto& distinct_count : _bin_distinct_counts) {
          _total_distinct_count += distinct_count;
        }
      }

template <typename T>
std::string MaxDiffFrHistogram<T>::name() const {
  return "MaxDiffFr";
}

struct ValueDistance {
  // The distance between the index-th and index+1-th element.
  uint32_t index;
  float distance;
};

bool sortDistance (ValueDistance el1, ValueDistance el2) {
  return (el1.distance > el2.distance);
}

bool sort_distances_per_index (ValueDistance el1, ValueDistance el2) {
  return (el1.index < el2.index);
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

      return std::make_shared<MaxDiffFrHistogram<T>>(std::move(bin_minima), 
          std::move(bin_maxima), 
          std::move(bin_heights), 
          std::move(bin_distinct_counts), total_count);
  }

  // Find the number of bins.
  auto bin_count = max_bin_count;
  if (static_cast<BinID>(value_distribution.size()) < max_bin_count) {
    bin_count = static_cast<BinID>(value_distribution.size());
  }

  std::vector<ValueDistance> distances(static_cast<int>(value_distribution.size() - 1));
  for (auto ind = uint32_t{0}; ind < value_distribution.size() - 1; ++ind) {
    struct ValueDistance val_dist;
    val_dist.index = ind;
    val_dist.distance = std::abs(value_distribution[ind].second - value_distribution[ind + 1].second);
    distances[ind] = val_dist;
  }

  std::sort(distances.begin(), distances.end(), sortDistance);
  distances.resize(bin_count - 1);
  std::sort(distances.begin(), distances.end(), sort_distances_per_index);

  std::vector<T> bin_minima(bin_count);
  std::vector<T> bin_maxima(bin_count);
  std::vector<HistogramCountType> bin_heights(bin_count, 0);
  std::vector<HistogramCountType> bin_distinct_counts(bin_count, 0);

  bin_minima[0] = value_distribution[0].first;
  bin_maxima[bin_count - 1] = value_distribution[value_distribution.size() - 1].first;

  for (auto bin_index = BinID{0}; bin_index < bin_count; ++bin_index) {
    if (bin_index > 0) {
      // Index of the next barrier position.
      const auto value_position = distances[bin_index - 1].index + 1;
      bin_minima[bin_index] = value_distribution[value_position].first;
    }
    if (bin_index < bin_count - 1) {
      const auto value_position = distances[bin_index].index;
      bin_maxima[bin_index] = value_distribution[value_position].first;
    }
  }

  auto bin_index = BinID{0};
  for (auto value_index = uint32_t{0}; value_index < value_distribution.size(); ++value_index) {
    const auto& [value, frequency] = value_distribution[value_index];

    while (value > bin_maxima[bin_index]) ++bin_index;

    bin_heights.at(bin_index) += frequency;
    ++bin_distinct_counts[bin_index];
  }

  return std::make_shared<MaxDiffFrHistogram<T>>(std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights), std::move(bin_distinct_counts), total_count);
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
  return _total_distinct_count;
}

template <typename T>
HistogramCountType MaxDiffFrHistogram<T>::bin_distinct_count(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");
  return _bin_distinct_counts[index];
}

EXPLICITLY_INSTANTIATE_DATA_TYPES(MaxDiffFrHistogram);

}  // namespace hyrise
