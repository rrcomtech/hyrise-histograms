#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tsl/robin_map.h"

#include "abstract_histogram.hpp"
#include "types.hpp"

namespace hyrise {

class Table;

template <typename T>
using ValueDistributionMap =
        tsl::robin_map<T, HistogramCountType, std::hash<T>, std::equal_to<T>,
                std::allocator<std::pair<T, HistogramCountType>>, std::is_same_v<std::decay_t<T>, pmr_string>>;

// TODO: Add an include guard. See https://en.wikipedia.org/wiki/Include_guard.
struct error_increase {
    uint32_t barrier_index;
    uint32_t error_increase;
};

struct error_decrease {
    uint32_t bin_index;
    uint32_t ideal_barrier_index;
    uint32_t error_decrease;
};

/**
 * Distinct-balanced histogram.
 * The first `bin_count_with_extra_value` contain each contain `distinct_count_per_bin + 1` distinct values,
 * all other bins contain `distinct_count_per_bin` distinct values.
 * There might be gaps between bins.
 */
template <typename T>
class GDYHistogram : public AbstractHistogram<T> {
 public:
  using AbstractHistogram<T>::AbstractHistogram;

  GDYHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                              std::vector<HistogramCountType>&& bin_heights,
                              const HistogramCountType distinct_count_per_bin, const BinID bin_count_with_extra_value,
                              const HistogramDomain<T>& domain = {});

  /**
   * Create an GDYHistogram for a column (spanning all Segments) of a Table
   * @param max_bin_count   Desired number of bins. Less might be created, but never more. Must not be zero.
   */
  static std::shared_ptr<GDYHistogram<T>> from_column(const Table& table, const ColumnID column_id,
                                                                     const BinID max_bin_count,
                                                                     const HistogramDomain<T>& domain = {});

  std::string name() const override;
  std::shared_ptr<AbstractHistogram<T>> clone() const override;
  HistogramCountType total_distinct_count() const override;
  HistogramCountType total_count() const override;

  /**
   * Returns the number of bins actually present in the histogram.
   * This number can be smaller than the number of bins requested when creating a histogram.
   * The number of bins is capped at the number of distinct values in the segment.
   * Otherwise, there would be empty bins without any benefit.
   */
  BinID bin_count() const override;

  const T& bin_minimum(const BinID index) const override;
  const T& bin_maximum(const BinID index) const override;
  HistogramCountType bin_height(const BinID index) const override;
  HistogramCountType bin_distinct_count(const BinID index) const override;

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

  // Number of distinct values per bin.
  HistogramCountType _distinct_count_per_bin;

  // The first bin_count_with_extra_value bins have an additional distinct value.
  BinID _bin_count_with_extra_value;

  // Aggregated counts over all bins, to avoid redundant computation
  HistogramCountType _total_count;
  HistogramCountType _total_distinct_count;

  // Returns the error increase, if a barrier would be removed.
  uint32_t error_increase_for_barrier_removal(uint32_t barrier_index);
  /*
   * Returns the error decrease for the ideal barrier insertion at the given index
   * and the index where the barrier would be inserted.
   * E.g.
   * error_decrease_for_barrier_insertion(1) = (5, 1)
   * 0 1 2 | 3 4 5 | 6 7 8 9
   *            |
   *         insert Barrier there to decrease error by 5
   */
  std::tuple<uint32_t, uint32_t> error_decrease_for_barrier_insertion(uint32_t segment_index);

  // Returns the error for a segment.
  uint32_t  error_for_segment(ValueDistributionMap<T>& value_distribution, uint32_t segment_index);

  // Populates error decreases and increases
  /*
  static std::tuple<std::vector<error_increase>, std::vector<error_decrease>> calculate_error_changes(
          std::vector<std::pair<T, HistogramCountType>> value_distribution,
          std::vector<T> bin_minima, std::vector<T> bin_maxima
  );*/
};

}  // namespace hyrise
