#include <limits>
#include <memory>
#include <string>

#include "base_test.hpp"

#include "statistics/statistics_objects/equi_height_histogram.hpp"
#include "statistics/statistics_objects/generic_histogram.hpp"
#include "utils/load_table.hpp"

namespace hyrise {

class EquiHeightHistogramTest : public BaseTest {
  void SetUp() override {
    _int_float4 = load_table("resources/test_data/tbl/int_float4.tbl");
    _float2 = load_table("resources/test_data/tbl/float2.tbl");
    _string2 = load_table("resources/test_data/tbl/string2.tbl");
  }

 protected:
  std::shared_ptr<Table> _int_float4;
  std::shared_ptr<Table> _float2;
  std::shared_ptr<Table> _string2;
};

TEST_F(EquiHeightHistogramTest, FromColumnString) {
    //StringHistogramDomain default_domain;
    const auto default_domain_histogram =
          EquiHeightHistogram<pmr_string>::from_column(*_string2, ColumnID{0}, 4u);

    ASSERT_EQ(default_domain_histogram->bin_count(), 4u);
    EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<pmr_string>("aa", "bla", 4, 4));
    EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<pmr_string>("bla", "uuu", 4, 4));
    EXPECT_EQ(default_domain_histogram->bin(BinID{2}), HistogramBin<pmr_string>("uuu", "yyy", 4, 4));
    EXPECT_EQ(default_domain_histogram->bin(BinID{2}), HistogramBin<pmr_string>("zzz", "zzz", 3, 3));
}

TEST_F(EquiHeightHistogramTest, FromColumnInt) {
    const auto hist = EquiHeightHistogram<int32_t>::from_column(*_int_float4, ColumnID{0}, 2u);

    ASSERT_EQ(hist->bin_count(), 2u);

    EXPECT_EQ(hist->bin_height(BinID{0}), HistogramCountType{4});
    EXPECT_EQ(hist->bin_height(BinID{1}), HistogramCountType{3});

    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<int32_t>(12, 12345, 4, 4));
    EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<int32_t>(123456, 123456, 3, 3));
}

TEST_F(EquiHeightHistogramTest, FromColumnFloat) {
    auto hist = EquiHeightHistogram<float>::from_column(*_int_float4, ColumnID{1}, 3u);

    ASSERT_EQ(hist->bin_count(), 3u);

    EXPECT_EQ(hist->bin_height(BinID{0}), HistogramCountType{3});
    EXPECT_EQ(hist->bin_height(BinID{1}), HistogramCountType{2});
    EXPECT_EQ(hist->bin_height(BinID{2}), HistogramCountType{2});
    
    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<float>(350.7f, 457.7f, 3, 3));
    EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<float>(458.7f, 700.0f, 2, 2));
    EXPECT_EQ(hist->bin(BinID{2}), HistogramBin<float>(800.0f, 900.0f, 2, 2));
}

}  // namespace hyrise
