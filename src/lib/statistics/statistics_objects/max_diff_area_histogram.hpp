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

struct ValueDistance {
  // The distance between the index-th and index+1-th element.
  uint32_t index;
  float distance;
};

template <typename T>
class MaxDiffAreaHistogram : public AbstractHistogram<T> {
 public:
  using AbstractHistogram<T>::AbstractHistogram;

  MaxDiffAreaHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                      std::vector<HistogramCountType>&& bin_height, std::vector<HistogramCountType>&& bin_distinct_counts, const HistogramCountType total_count,
                      const HistogramCountType total_distinct_count, const HistogramDomain<T>& domain = {});

  static std::shared_ptr<MaxDiffAreaHistogram<T>> from_column(const Table& table, const ColumnID column_id,
                                                             const BinID max_bin_count,
                                                             const HistogramDomain<T>& domain = {});

  std::string name() const override;
  std::shared_ptr<AbstractHistogram<T>> clone() const override;
  HistogramCountType total_distinct_count() const override;
  HistogramCountType total_count() const override;

  BinID bin_count() const override;

  const T& bin_minimum(BinID index) const override;
  const T& bin_maximum(BinID index) const override;
  HistogramCountType bin_height(const BinID index) const override;
  HistogramCountType bin_distinct_count(const BinID index) const override;

  bool operator==(const MaxDiffAreaHistogram<T>& rhs) const;

  static bool sortDistance (ValueDistance el1, ValueDistance el2) {
    return (el1.distance > el2.distance);
  }

  static bool sort_distances_per_index (ValueDistance el1, ValueDistance el2) {
    return (el1.index < el2.index);
  }

 protected:
  BinID _bin_for_value(const T& value) const override;

  BinID _next_bin_for_value(const T& value) const override;

 private:
  /**
   * We use multiple vectors rather than a vector of structs for ease-of-use with STL library functions.
   */
  // Min values on a per-bin basis.
  std::vector<T> _bin_minima;

  // Max values on a per-bin basis.
  std::vector<T> _bin_maxima;

  // Number of values on a per-bin basis.
  std::vector<HistogramCountType> _bin_heights;

  // Number of distinct values on a per-bin basis.
  std::vector<HistogramCountType> _bin_distinct_counts;

  // Aggregated counts over all bins, to avoid redundant computation
  HistogramCountType _total_count;
  HistogramCountType _total_distinct_count = 0;
};

// For gtest
template <typename T>
std::ostream& operator<<(std::ostream& stream, const MaxDiffAreaHistogram<T>& histogram) {
  stream << histogram.description() << std::endl;
  return stream;
}

// EXPLICITLY_DECLARE_DATA_TYPES(MaxDiffAreaHistogram);

}  // namespace hyrise