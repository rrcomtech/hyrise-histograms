#include <limits>
#include <memory>
#include <string>

#include "base_test.hpp"
#include "statistics/statistics_objects/equal_distinct_count_histogram.hpp"
#include "statistics/statistics_objects/generic_histogram.hpp"
#include "utils/load_table.hpp"

namespace hyrise {

class EqualDistinctCountHistogramTest : public BaseTest {
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

TEST_F(EqualDistinctCountHistogramTest, Name) {
  const auto histogram = EqualDistinctCountHistogram<int32_t>{{1}, {100}, {50}, 10, 0};
  EXPECT_EQ(histogram.name(), "EqualDistinctCount");
}

TEST_F(EqualDistinctCountHistogramTest, FromColumnString) {
  const auto default_domain = StringHistogramDomain{};
  const auto default_domain_histogram =
      EqualDistinctCountHistogram<pmr_string>::from_column(*_string2, ColumnID{0}, 4, default_domain);

  ASSERT_EQ(default_domain_histogram->bin_count(), 4);
  EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<pmr_string>("aa", "birne", 3, 3));
  EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<pmr_string>("bla", "ttt", 4, 3));
  EXPECT_EQ(default_domain_histogram->bin(BinID{2}), HistogramBin<pmr_string>("uuu", "xxx", 4, 3));

  const auto reduced_histogram = StringHistogramDomain{'a', 'c', 9};
  const auto reduced_domain_histogram =
      EqualDistinctCountHistogram<pmr_string>::from_column(*_string2, ColumnID{0}, 4, reduced_histogram);

  ASSERT_EQ(default_domain_histogram->bin_count(), 4);
  EXPECT_EQ(default_domain_histogram->bin(BinID{0}), HistogramBin<pmr_string>("aa", "birne", 3, 3));
  EXPECT_EQ(default_domain_histogram->bin(BinID{1}), HistogramBin<pmr_string>("bla", "ttt", 4, 3));
  EXPECT_EQ(default_domain_histogram->bin(BinID{2}), HistogramBin<pmr_string>("uuu", "xxx", 4, 3));
}

TEST_F(EqualDistinctCountHistogramTest, FromColumnInt) {
  const auto hist = EqualDistinctCountHistogram<int32_t>::from_column(*_int_float4, ColumnID{0}, 2);

  ASSERT_EQ(hist->bin_count(), 2);
  EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<int32_t>(12, 123, 2, 2));
  EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<int32_t>(12345, 123456, 5, 2));
}

TEST_F(EqualDistinctCountHistogramTest, FromColumnFloat) {
  const auto hist = EqualDistinctCountHistogram<float>::from_column(*_float2, ColumnID{0}, 3);

  ASSERT_EQ(hist->bin_count(), 3);
  EXPECT_EQ(hist->bin(BinID{0}), HistogramBin<float>(0.5f, 2.2f, 4, 4));
  EXPECT_EQ(hist->bin(BinID{1}), HistogramBin<float>(2.5f, 3.3f, 6, 3));
  EXPECT_EQ(hist->bin(BinID{2}), HistogramBin<float>(3.6f, 6.1f, 4, 3));
}

}  // namespace hyrise
