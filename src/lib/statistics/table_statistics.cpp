#include "table_statistics.hpp"

#include <algorithm>
#include <numeric>
#include <thread>
#include <cstdlib>

#include "attribute_statistics.hpp"
#include "hyrise.hpp"
#include "resolve_type.hpp"
#include "scheduler/job_task.hpp"
#include "statistics/statistics_objects/abstract_histogram.hpp"
#include "statistics/statistics_objects/equi_height_histogram.hpp"
#include "statistics/statistics_objects/equi_width_histogram.hpp"
#include "statistics/statistics_objects/max_diff_fr_histogram.hpp"
#include "storage/table.hpp"
#include "utils/assert.hpp"

namespace hyrise {

std::shared_ptr<TableStatistics> TableStatistics::from_table(const Table& table) {
  std::vector<std::shared_ptr<BaseAttributeStatistics>> column_statistics(table.column_count());

  /**
   * Determine bin count, within mostly arbitrarily chosen bounds: 5 (for tables with <=2k rows) up to 100 bins
   * (for tables with >= 200m rows) are created.
   */
  const auto histogram_bin_count = std::min<size_t>(100, std::max<size_t>(5, table.row_count() / 2'000));

  /**
   * We highly recommend setting up a multithreaded scheduler before the following procedure is executed to parallelly
   * create statistics objects for the table's columns.
   */

  auto jobs = std::vector<std::shared_ptr<AbstractTask>>{};
  jobs.reserve(table.column_count());

  for (auto column_id = ColumnID{0}; column_id < table.column_count(); ++column_id) {
    const auto generate_column_statistics = [&, column_id]() {
      const auto column_data_type = table.column_data_type(column_id);
      resolve_data_type(column_data_type, [&](auto type) {
        using ColumnDataType = typename decltype(type)::type;

        const auto output_column_statistics = std::make_shared<AttributeStatistics<ColumnDataType>>();
        std::shared_ptr<AbstractHistogram<ColumnDataType>> histogram;

        // The type of Histogram, that should be used. Set it in the Environment Variable by:
        // export HISTOGRAM=<type>
        const auto HISTORAM_TYPE = std::getenv("HISTOGRAM");

        if (HISTORAM_TYPE) {
            if (strcmp(HISTORAM_TYPE, "EquiHeightHistogram") == 0) {
                PerformanceWarning("EquiHeightHistogram is used!");
                histogram = EquiHeightHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
            } else if (strcmp(HISTORAM_TYPE, "EquiWidthHistogram") == 0) {
              // EquiWidth cannot handle strings well. Therefore, EquiHeight will be used in that case.
              if constexpr (std::is_same_v<ColumnDataType, pmr_string>) {
                PerformanceWarning("Fallback EquiHeightHistogram is used!");
                histogram = EquiHeightHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
              } else {
                PerformanceWarning("EquiWidthHistogram is used!");
                histogram = EquiWidthHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
              }
            } else if (strcmp(HISTORAM_TYPE, "EqualDistinctCountHistogram") == 0) {
                PerformanceWarning("EqualDistinctCountHistogram is used!");
                histogram = EqualDistinctCountHistogram<ColumnDataType>::from_column(table, column_id,
                                                                                     histogram_bin_count);
            } else if (strcmp(HISTORAM_TYPE, "MaxDiffHistogram") == 0) {
                if constexpr (std::is_same_v<ColumnDataType, pmr_string>) {
                    PerformanceWarning("Fallback EquiHeightHistogram is used!");
                    histogram = EquiHeightHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
                } else {
                    PerformanceWarning("MaxDiffHistogram (Frequency-based) is used!");
                    histogram = MaxDiffFrHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
                }
            } else {
              Fail("Unknown Histogram specified!");
            }
        } else {
          PerformanceWarning("EqualDistinctCountHistogram is used!");
          histogram = EqualDistinctCountHistogram<ColumnDataType>::from_column(table, column_id, histogram_bin_count);
        }

        if (histogram) {
          output_column_statistics->set_statistics_object(histogram);

          // Use the insight that the histogram will only contain non-null values to generate the NullValueRatio
          // property
          const auto null_value_ratio =
              table.row_count() == 0
                  ? 0.0f
                  : 1.0f - (static_cast<float>(histogram->total_count()) / static_cast<float>(table.row_count()));
          output_column_statistics->set_statistics_object(std::make_shared<NullValueRatioStatistics>(null_value_ratio));
        } else {
          // Failure to generate a histogram currently only stems from all-null segments.
          // TODO(anybody) this is a slippery assumption. But the alternative would be a full segment scan...
          output_column_statistics->set_statistics_object(std::make_shared<NullValueRatioStatistics>(1.0f));
        }

        column_statistics[column_id] = output_column_statistics;
      });
    };
    jobs.emplace_back(std::make_shared<JobTask>(generate_column_statistics));
  }
  Hyrise::get().scheduler()->schedule_and_wait_for_tasks(jobs);

  return std::make_shared<TableStatistics>(std::move(column_statistics), table.row_count());
}

TableStatistics::TableStatistics(std::vector<std::shared_ptr<BaseAttributeStatistics>>&& init_column_statistics,
                                 const Cardinality init_row_count)
    : column_statistics(std::move(init_column_statistics)), row_count(init_row_count) {}

DataType TableStatistics::column_data_type(const ColumnID column_id) const {
  DebugAssert(column_id < column_statistics.size(), "ColumnID out of bounds");
  return column_statistics[column_id]->data_type;
}

std::ostream& operator<<(std::ostream& stream, const TableStatistics& table_statistics) {
  stream << "TableStatistics {" << std::endl;
  stream << "  RowCount: " << table_statistics.row_count << "; " << std::endl;

  for (const auto& column_statistics : table_statistics.column_statistics) {
    resolve_data_type(column_statistics->data_type, [&](const auto data_type_t) {
      using ColumnDataType = typename decltype(data_type_t)::type;
      stream << *std::dynamic_pointer_cast<AttributeStatistics<ColumnDataType>>(column_statistics) << std::endl;
    });
  }

  stream << "}" << std::endl;

  return stream;
}

}  // namespace hyrise
