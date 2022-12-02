#include <limits>
#include <memory>
#include <string>

#include "base_test.hpp"

#include "statistics/statistics_objects/equi_width_histogram.hpp"
#include "statistics/statistics_objects/generic_histogram.hpp"
#include "utils/load_table.hpp"

namespace hyrise {

class EquiWidthHistogramTest : public BaseTest {
  void SetUp() override {
    _int_float4 = load_table("resources/test_data/tbl/int_float4.tbl");
    _float2 = load_table("resources/test_data/tbl/float2.tbl");
    _int_highly_frequent = load_table("resources/test_data/tbl/int_highly_frequent.tbl");
  }

 protected:
  std::shared_ptr<Table> _int_float4;
  std::shared_ptr<Table> _float2;
  std::shared_ptr<Table> _int_highly_frequent;
};

TEST_F(EquiWidthHistogramTest, FromHighlyFrequentInt) {
    //StringHistogramDomain default_domain;
    const auto default_domain_histogram =
            EquiWidthHistogram<int32_t>::from_column(*_int_highly_frequent, ColumnID{0}, 3u);

    ASSERT_EQ(default_domain_histogram->bin_count(), 3u);
    EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<int32_t>(1, 1, 1, 1));
    EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<int32_t>(2, 2, 1, 1));
    EXPECT_EQ(default_domain_histogram->bin(BinID{2}), HistogramBin<int32_t>(3, 3, 16, 1));

}

TEST_F(EquiWidthHistogramTest, FromColumnInt) {
    const auto hist = EquiWidthHistogram<int32_t>::from_column(*_int_float4, ColumnID{0}, 2u);

    ASSERT_EQ(hist->bin_count(), 2u);
    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<int32_t>(12, 12345, 4, 3));
    EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<int32_t>(123456, 123456, 3, 1));
}

TEST_F(EquiWidthHistogramTest, FromColumnFloat) {
    auto hist = EquiWidthHistogram<float>::from_column(*_int_float4, ColumnID{1}, 3u);

    ASSERT_EQ(hist->bin_count(), 3u);
    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<float>(350.7f, 457.7f, 3, 3));
    EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<float>(458.7f, 800.0f, 3, 3));
    EXPECT_EQ(hist->bin(BinID{2}), HistogramBin<float>(900.0f, 900.0f, 1, 1));
}

}  // namespace hyrise
