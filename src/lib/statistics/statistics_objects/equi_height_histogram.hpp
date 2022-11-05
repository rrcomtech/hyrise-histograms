#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tsl/robin_map.h"

#include "abstract_histogram.hpp"
#include "types.hpp"

namespace hyrise {

template <typename T>
class EquiHeigthHistogram : public AbstractHistogram<T> {
 public:
  using AbstractHistogram<T>::AbstractHistogram;

  EquiHeigthHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                   HistogramCountType&& bin_height, std::vector<HistogramCountType>&& bin_distinct_counts,
                   HistogramCountType&& total_count, const HistogramDomain<T>& domain = {});

  // Convenience builder for a EquiHeigthHistogram with a single bin
  static std::shared_ptr<EquiHeigthHistogram<T>> with_single_bin(const T& min, const T& max,
                                                              const HistogramCountType& height,
                                                              const HistogramCountType& distinct_count,
                                                              const HistogramDomain<T>& domain = {});

  static std::shared_ptr<EquiHeigthHistogram<T>> from_column(const Table& table, const ColumnID column_id,
                                                                     const BinID max_bin_count,
                                                                     const HistogramDomain<T>& domain = {});

  std::string name() const override;
  HistogramCountType total_count() const override;

  BinID bin_count() const override;

  const T& bin_minimum(const BinID index) const override;
  const T& bin_maximum(const BinID index) const override;
  HistogramCountType bin_height(const BinID index) const override;

  bool operator==(const EquiHeigthHistogram<T>& rhs) const;

 protected:
  using ValueDistributionMap =
    tsl::robin_map<T, HistogramCountType, std::hash<T>, std::equal_to<T>,
                   std::allocator<std::pair<T, HistogramCountType>>, std::is_same_v<std::decay_t<T>, pmr_string>>;

  BinID _bin_for_value(const T& value) const override;

  BinID _next_bin_for_value(const T& value) const override;

  void _add_segment_to_value_distribution(const AbstractSegment& segment, ValueDistributionMap& value_distribution,
                                       const HistogramDomain<T>& domain);

  std::vector<std::pair<T, HistogramCountType>> _value_distribution_from_column(const Table& table,
                                                                             const ColumnID column_id,
                                                                             const HistogramDomain<T>& domain);

 private:
  /**
   * We use multiple vectors rather than a vector of structs for ease-of-use with STL library functions.
   */
  // Min values on a per-bin basis.
  std::vector<T> _bin_minima;

  // Max values on a per-bin basis.
  std::vector<T> _bin_maxima;

  // Number of values on a per-bin basis.
  HistogramCountType _bin_height;

  // Number of distinct values on a per-bin basis.
  std::vector<HistogramCountType> _bin_distinct_counts;

  // Aggregated counts over all bins, to avoid redundant computation
  HistogramCountType _total_count;
};

// For gtest
template <typename T>
std::ostream& operator<<(std::ostream& stream, const EquiHeigthHistogram<T>& histogram) {
  stream << histogram.description() << std::endl;
  return stream;
}

EXPLICITLY_DECLARE_DATA_TYPES(EquiHeigthHistogram);

}  // namespace hyrise
