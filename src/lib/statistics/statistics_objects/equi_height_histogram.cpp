#include "equi_height_histogram.hpp"

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

namespace hyrise {

    template<typename T>
    EquiHeightHistogram<T>::EquiHeightHistogram(std::vector<T> &&bin_minima, std::vector<T> &&bin_maxima,
                                                HistogramCountType &&bin_height, HistogramCountType &&total_count,
                                                const HistogramDomain<T> &domain) :
            _bin_minima(std::move(bin_minima)),
            _bin_maxima(std::move(bin_maxima)),
            _bin_heights(std::move(bin_height)),
            _total_count(total_count) {}

    template<typename T>
    std::string EquiHeightHistogram<T>::name() const {
        return "Equi Heigth";
    }

    template<typename T>
    std::shared_ptr<EquiHeightHistogram<T>>
    EquiHeightHistogram<T>::from_column(const Table &table, const ColumnID column_id,
                                        const BinID max_bin_count,
                                        const HistogramDomain<T> &domain) {
        Assert(max_bin_count > 0, "max_bin_count must be greater than zero ");

        const auto value_distribution = value_distribution_from_column(table, column_id, domain);

        if (value_distribution.empty()) {
            return nullptr;
        }

        const auto adder = [](const std::pair<T, HistogramCountType> &a, const std::pair<T, HistogramCountType> &b) {
            return a.second + b.second;
        };

        const auto total_count = std::accumulate(value_distribution.begin(), value_distribution.end(), 0, adder);
        const auto bin_count = total_count ? total_count < max_bin_count : max_bin_count;
        const auto values_per_bin = std::ceil(total_count / bin_count);

        std::vector<T> bin_minima(bin_count);
        std::vector<T> bin_maxima(bin_count);
        std::vector<T> bin_heights(bin_count);

        const auto sorter = [](const std::pair<T, HistogramCountType> &a, const std::pair<T, HistogramCountType> &b) {
            return a.first < b.first;
        };

        std::sort(value_distribution.begin(), value_distribution.end(), sorter);

        auto value_distribution_index = 0;
        for (auto bin_id = BinID{0}; bin_id < bin_count; ++bin_id) {
            bin_minima[bin_id] = value_distribution[value_distribution_index].first;
            auto space_left_in_bin = values_per_bin;
            auto values_left_for_value = value_distribution[value_distribution_index].second;

            // Do until every bucket is reached and
            while (space_left_in_bin > 0 && value_distribution_index < value_distribution.size()) {

                if (space_left_in_bin == values_left_for_value) {
                    // The amount of values exactly fills the bucket.
                    bin_maxima[bin_id] = value_distribution[value_distribution_index].first;
                    bin_heights[bin_id] = values_per_bin;
                    space_left_in_bin -= values_left_for_value;
                    ++value_distribution_index;
                } else {
                    if (space_left_in_bin > values_left_for_value) {
                        // Fewer values than space in bin.
                        space_left_in_bin -= values_left_for_value;
                        ++value_distribution_index;
                        if (value_distribution.size() == 0) {
                            // If bin is last bin and can't be fully filled, because no more values present.
                            bin_maxima[bin_id] = value_distribution[value_distribution_index].first;
                            bin_heights[bin_id] = values_per_bin - space_left_in_bin;
                        }
                    } else {
                        // There are more values than space in bin.
                        values_left_for_value -= space_left_in_bin;
                        space_left_in_bin = 0;
                        ++value_distribution_index;
                        bin_maxima[bin_id] = value_distribution[value_distribution_index].first;
                        bin_heights[bin_id] = values_per_bin;
                    }
                }
            }
        }

        return std::make_shared<EquiHeightHistogram<T>>(
                std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
                static_cast<HistogramCountType>(bin_count), values_per_bin);
    }

    template<typename T>
    void EquiHeightHistogram<T>::_add_segment_to_value_distribution(const AbstractSegment &segment,
                                                                    ValueDistributionMap &value_distribution,
                                                                    const HistogramDomain<T> &domain) {
        segment_iterate<T>(segment, [&](const auto &iterator_value) {
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

    template<typename T>
    std::vector<std::pair<T, HistogramCountType>>
    EquiHeightHistogram<T>::_value_distribution_from_column(const Table &table,
                                                            const ColumnID column_id,
                                                            const HistogramDomain<T> &domain) {
        auto value_distribution_map = ValueDistributionMap{};
        const auto chunk_count = table.chunk_count();
        for (auto chunk_id = ChunkID{0}; chunk_id < chunk_count; ++chunk_id) {
            const auto chunk = table.get_chunk(chunk_id);
            if (!chunk) {
                continue;
            }

            add_segment_to_value_distribution<T>(*chunk->get_segment(column_id), value_distribution_map, domain);
        }

        auto value_distribution =
                std::vector<std::pair<T, HistogramCountType>>{value_distribution_map.begin(),
                                                              value_distribution_map.end()};
        value_distribution_map.clear();  // Maps can be large and sorting slow. Free space early.
        boost::sort::pdqsort(value_distribution.begin(), value_distribution.end(),
                             [&](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

        return value_distribution;
    }

    template<typename T>
    HistogramCountType EquiHeightHistogram<T>::total_count() const {
        return _total_count;
    }

    template<typename T>
    BinID EquiHeightHistogram<T>::bin_count() const {
        return _bin_maxima.size();
    }

    template<typename T>
    const T& EquiHeightHistogram<T>::bin_minimum(BinID index) const {
        return _bin_minima[index];
    }

    template<typename T>
    const T& EquiHeightHistogram<T>::bin_maximum(BinID index) const {
        return _bin_maxima[index];
    }
/*
    template<typename T>
    HistogramCountType EquiHeightHistogram<T>::bin_height(const BinID index) const {
        return _bin_heights[index];
    }
*/
    template<typename T>
    BinID EquiHeightHistogram<T>::_bin_for_value(const T& value) const {
        for (int index = 0; index < this->bin_count() - 1; ++index) {
            if (bin_minimum(index + 1) > value) {
                return index;
            }
        }
        // Value must be in last bin.
        return this->bin_count() - 1;
    }

    template<typename T>
    BinID EquiHeightHistogram<T>::_next_bin_for_value(const T& value) const {
        const auto bin_for_value = this->_bin_for_value(value);
        // If value is in last bin, return first bin.
        if (bin_for_value == bin_count() - 1) {
            return 0;
        }
        return bin_for_value + 1;
    }

}