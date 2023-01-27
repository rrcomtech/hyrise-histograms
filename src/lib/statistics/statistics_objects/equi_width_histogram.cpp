#include "equi_width_histogram.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <typeinfo>

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
                   std::allocator<std::pair<T, HistogramCountType>>>;

template <typename T>
void add_segment_to_value_distribution(const AbstractSegment& segment, ValueDistributionMap<T>& value_distribution,
                                       const HistogramDomain<T>& domain) {
  segment_iterate<T>(segment, [&](const auto& iterator_value) {
    if (iterator_value.is_null()) {
      return;
    }
      ++value_distribution[iterator_value.value()];
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
EquiWidthHistogram<T>::EquiWidthHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                                            std::vector<HistogramCountType>&& bin_height, std::vector<HistogramCountType>&& bin_distinct_counts,
                                            const HistogramCountType total_count, const HistogramDomain<T>& domain)
    : AbstractHistogram<T>(domain),
    _bin_minima(std::move(bin_minima)),
      _bin_maxima(std::move(bin_maxima)),
      _bin_heights(std::move(bin_height)),
      _bin_distinct_counts(std::move(bin_distinct_counts)),
      _total_count{total_count} {
        _total_distinct_count = std::accumulate(_bin_distinct_counts.cbegin(), _bin_distinct_counts.cend(), HistogramCountType{0});
      }

template <typename T>
std::string EquiWidthHistogram<T>::name() const {
  return "EquiWidth";
}


template <typename T>
std::shared_ptr<EquiWidthHistogram<T>> EquiWidthHistogram<T>::from_column(const Table& table,
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

  const auto max_value = value_distribution[value_distribution.size() - 1].first;
  const auto min_value = value_distribution[0].first;
  auto bin_width = (static_cast<double>(max_value) - static_cast<double>(min_value)) / static_cast<double>(bin_count);

  // Initialize the resulting data structures.
  std::vector<T> bin_minima(bin_count, max_value);
  std::vector<T> bin_maxima(bin_count, min_value);

  std::vector<bool> minimum_set(bin_count, false);
  std::vector<bool> maximum_set(bin_count, false);

  std::vector<HistogramCountType> bin_heights(bin_count, 0);
  std::vector<HistogramCountType> bin_distinct_counts(bin_count, 0);

  if (bin_width == 0) {
    Assert(min_value == max_value, "Invalid bin_width without reason!");
    ++bin_width;  
  }

  BinID bin_index;
  const auto value_distribution_size = value_distribution.size();
  for (auto value_index = size_t{0}; value_index < value_distribution_size; ++value_index) {
    const auto value = value_distribution[value_index].first;
    const auto value_from_min_value = value - min_value;

    if (value == max_value) {
      bin_index = bin_count - 1;
    } else {
      bin_index = static_cast<uint32_t>(std::floor(static_cast<double>(value_from_min_value) / bin_width));
    }


    Assert(bin_index < bin_minima.size(), "bin_index out of range!");

    bin_heights[bin_index] += value_distribution[value_index].second;
    ++bin_distinct_counts[bin_index];

    if (bin_minima[bin_index] > value) {
      bin_minima[bin_index] = value;
      minimum_set[bin_index] = true;
    }
    if (bin_maxima[bin_index] < value) {
      bin_maxima[bin_index] = value;
      maximum_set[bin_index] = true;
    }
    
    if (bin_maxima[bin_index] < bin_minima[bin_index]) {
      std::cout << "\n\n\n\n\nSetting wrond bin!\n" << std::endl;
      std::cout << "value: " << value << "\n min_value: " << min_value << "\nmax_value: " << max_value << "\nbin_width: " << bin_width << std::endl;
      std::cout << "value - min_value: " << value_from_min_value << "\nbin_index: " << bin_index << "\nbin_count: " << bin_count << std::endl;
      std::cout << "bin_min: " << bin_minima[bin_index] << " bin_max: " << bin_maxima[bin_index] << std::endl;
      Fail("Invalid Bin in EquiWidth Construction!");
    }

  }

  for (auto index = size_t{0}; index < maximum_set.size(); ++index) {
    if (!minimum_set[index]) {
      bin_minima[index] = bin_maxima.at(index - 1);
    }
    if (!maximum_set[index]) {
      bin_maxima[index] = bin_minima.at(index);
    }
  }

  return std::make_shared<EquiWidthHistogram<T>>(std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights), std::move(bin_distinct_counts), total_count);
}

template <typename T>
HistogramCountType EquiWidthHistogram<T>::total_count() const {
  return _total_count;
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> EquiWidthHistogram<T>::clone() const {
  // The new histogram needs a copy of the data
  auto bin_minima_copy = _bin_minima;
  auto bin_maxima_copy = _bin_maxima;
  auto bin_heights_copy = _bin_heights;
  auto bin_distinct_counts_copy = _bin_distinct_counts;

  return std::make_shared<EquiWidthHistogram<T>>(std::move(bin_minima_copy), std::move(bin_maxima_copy),
                                                  std::move(bin_heights_copy), std::move(bin_distinct_counts_copy), _total_count);
}

template <typename T>
BinID EquiWidthHistogram<T>::bin_count() const {
  return _bin_maxima.size();
}

template <typename T>
const T& EquiWidthHistogram<T>::bin_minimum(BinID index) const {
  return _bin_minima.at(index);
}

template <typename T>
const T& EquiWidthHistogram<T>::bin_maximum(BinID index) const {
  // auto max = _bin_maxima.at(index);
  // if (max == 0 && bin_minimum(index) > 0) {
  //   std::cout << "IT HAPPPPPPPPPPPENEDEDED ##################" << std::endl;
  // }

  return _bin_maxima.at(index);
}

template <typename T>
HistogramCountType EquiWidthHistogram<T>::bin_height(const BinID index) const {
  return _bin_heights[index];
}

template <typename T>
BinID EquiWidthHistogram<T>::_bin_for_value(const T& value) const {
  for (auto index = size_t{0}; index < this->bin_count() - 1; ++index) {
    if (bin_minimum(index + 1) > value) {
      return index;
    }
  }
  // Value must be in last bin.
  return this->bin_count() - 1;
}

template <typename T>
BinID EquiWidthHistogram<T>::_next_bin_for_value(const T& value) const {
  const auto bin_for_value = this->_bin_for_value(value);
  // If value is in last bin, return first bin.
  if (bin_for_value == bin_count() - 1) {
    return 0;
  }
  return bin_for_value + 1;
}

template <typename T>
HistogramCountType EquiWidthHistogram<T>::total_distinct_count() const {
  return _total_distinct_count;
}

template <typename T>
HistogramCountType EquiWidthHistogram<T>::bin_distinct_count(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");
  return _bin_distinct_counts[index];
}

template class EquiWidthHistogram<int>;
template class EquiWidthHistogram<long>;
template class EquiWidthHistogram<float>;
template class EquiWidthHistogram<double>;
 
}  // namespace hyrise
