#include <limits>
#include <memory>
#include <string>

#include "base_test.hpp"

#include "statistics/statistics_objects/max_diff_fr_histogram.hpp"
#include "statistics/statistics_objects/generic_histogram.hpp"
#include "utils/load_table.hpp"

namespace hyrise {

class MaxDiffFrHistogramTest : public BaseTest {
  void SetUp() override {
    _int_float4 = load_table("resources/test_data/tbl/int_float4.tbl");
    _float2 = load_table("resources/test_data/tbl/float2.tbl");
    _string2 = load_table("resources/test_data/tbl/string2.tbl");
    _int_highly_frequent = load_table("resources/test_data/tbl/int_highly_frequent.tbl");
  }

 protected:
  std::shared_ptr<Table> _int_float4;
  std::shared_ptr<Table> _float2;
  std::shared_ptr<Table> _string2;
  std::shared_ptr<Table> _int_highly_frequent;
};

TEST_F(MaxDiffFrHistogramTest, FromHighlyFrequentInt) {
    //StringHistogramDomain default_domain;
    const auto default_domain_histogram =
          MaxDiffFrHistogram<int32_t>::from_column(*_int_highly_frequent, ColumnID{0}, 3u);
      
    ASSERT_EQ(default_domain_histogram->bin_count(), 2u);
    EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<int32_t>(1, 2, 2, 2));
    EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<int32_t>(3, 3, 16, 1));

}

TEST_F(MaxDiffFrHistogramTest, FromColumnString) {
    //StringHistogramDomain default_domain;
    const auto default_domain_histogram =
          MaxDiffFrHistogram<pmr_string>::from_column(*_string2, ColumnID{0}, 4u);

    ASSERT_EQ(default_domain_histogram->bin_count(), 2u);
    EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<pmr_string>("aa", "yyy", 12, 10));
    EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<pmr_string>("zzz", "zzz", 3, 1));
}

TEST_F(MaxDiffFrHistogramTest, FromColumnInt) {
    const auto hist = MaxDiffFrHistogram<int32_t>::from_column(*_int_float4, ColumnID{0}, 2u);

    ASSERT_EQ(hist->bin_count(), 1u);
    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<int32_t>(12, 123456, 7, 4));
}

TEST_F(MaxDiffFrHistogramTest, FromColumnFloat) {
    auto hist = MaxDiffFrHistogram<float>::from_column(*_int_float4, ColumnID{1}, 3u);

    ASSERT_EQ(hist->bin_count(), 1u);
    EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<float>(350.7f, 900.0f, 7, 7));
}

}  // namespace hyrise
