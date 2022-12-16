#include "benchmark_sql_executor.hpp"

#include "benchmark_runner.hpp"
#include "hyrise.hpp"
#include "operators/table_scan.hpp"
#include "operators/aggregate_hash.hpp"
#include "operators/join_hash.hpp"
#include "sql/sql_pipeline_builder.hpp"
#include "utils/check_table_equal.hpp"
#include "utils/timer.hpp"
#include "visualization/lqp_visualizer.hpp"
#include "visualization/pqp_visualizer.hpp"
#include "operators/pqp_utils.hpp"

namespace hyrise {
BenchmarkSQLExecutor::BenchmarkSQLExecutor(const std::shared_ptr<SQLiteWrapper>& sqlite_wrapper,
                                           const std::optional<std::string>& visualize_prefix)
    : _sqlite_connection(sqlite_wrapper ? std::optional<SQLiteWrapper::Connection>{sqlite_wrapper->new_connection()}
                                        : std::optional<SQLiteWrapper::Connection>{}),
      _visualize_prefix(visualize_prefix) {
  if (_sqlite_connection) {
    _sqlite_connection->raw_execute_query("BEGIN TRANSACTION");
    _sqlite_transaction_open = true;
  }
}

std::pair<SQLPipelineStatus, std::shared_ptr<const Table>> BenchmarkSQLExecutor::execute(
    const std::string& sql, const std::shared_ptr<const Table>& expected_result_table) {
  auto pipeline_builder = SQLPipelineBuilder{sql};
  if (transaction_context) {
    pipeline_builder.with_transaction_context(transaction_context);
  }

  auto pipeline = pipeline_builder.create_pipeline();

  const auto [pipeline_status, result_table] = pipeline.get_result_table();

  // #####################
  // ##################### Martins Code
  // #####################

  auto cardinality_estimator = CardinalityEstimator{};
  cardinality_estimator.guarantee_bottom_up_construction();

  const auto& lqp_plans = pipeline.get_optimized_logical_plans();
  // if (lqp_plans.size() != 1) {
  //   std::cerr << "!!!!!!\n!!!!!!\n!!!!!!\n";
  //   std::cerr << "!!!!!! Caution: analyzing a pipeline with more than one statement.\n";
  //   std::cerr << "!!!!!!\n!!!!!!\n!!!!!!\n" << std::endl;
  // }

  // Prints the logical query plans. Helps for understanding what's actually happening in the query.
  // We also store visualizations of these plans here: https://hyrise-ci.epic-hpi.de/job/hyrise/job/hyrise/job/master/lastStableBuild/artifact/query_plans/tpch/
  // for (const auto& plan : lqp_plans) {
  //   std::cout << "\n\n\n\n\n###\n###\n###\n" << *plan << "###\n###\n###\n" << std::endl;
  // }

  const auto& pqp_plans = pipeline.get_physical_plans();
  Assert(pqp_plans.size() == lqp_plans.size(), "NOTTHESAME");
  for (const auto& plan_root : pqp_plans) {
    if (plan_root->type() == OperatorType::CreateView || plan_root->type() == OperatorType::DropView) {
      // Skip some case. You might find more cases that need to be discarded.
      std::cout << "Skipping create/drop view operators." << std::endl;
      continue;
    }

    auto visitor = [&](const auto& node) {
      if (node->type () == OperatorType::TableScan) {
        if (node->left_input()->type() == OperatorType::GetTable || node->left_input()->type() == OperatorType::Validate) {
          const auto table_scan_op = std::dynamic_pointer_cast<const TableScan>(node);
          const auto& table_scan_performance_data = table_scan_op->performance_data;
          const auto& get_table_performance_data = node->left_input()->performance_data;

          Hyrise::get().cardinality_statistics += Hyrise::get().current_benchmark;
          Hyrise::get().cardinality_statistics += std::to_string(Hyrise::get().current_operator);
          Hyrise::get().cardinality_statistics += ",TableScan,";
          Hyrise::get().cardinality_statistics += std::to_string(get_table_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(table_scan_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node)) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(table_scan_op->lqp_node)) + "\n";

          // std::cout << "Information for " << *node->lqp_node << std::endl;
          // const auto table_scan_op = std::dynamic_pointer_cast<const TableScan>(node);
          // const auto& table_scan_performance_data = table_scan_op->performance_data;
          // const auto& get_table_performance_data = node->left_input()->performance_data;
          // std::cout << "Input size to table scan: " << get_table_performance_data->output_row_count << std::endl;
          // std::cout << "Result size after table scan: " << table_scan_performance_data->output_row_count << std::endl;
          // // std::cout << "Actual selectivity of table scan: " << static_cast<float>(get_table_performance_data->output_row_count) / static_cast<float>(table_scan_performance_data->output_row_count) << std::endl;

          // std::cout << "Input size to table scan: " << cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node) << std::endl;
          // std::cout << "Result size after table scan: " << cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node) << std::endl;
        }
      } else if (node->type () == OperatorType::Aggregate) {
        if (node->left_input()->type() == OperatorType::GetTable || node->left_input()->type() == OperatorType::Validate) {
          const auto aggregate_op = std::dynamic_pointer_cast<const AggregateHash>(node);
          const auto& aggregate_performance_data = aggregate_op->performance_data;
          const auto& get_table_performance_data = node->left_input()->performance_data;

          Hyrise::get().cardinality_statistics += Hyrise::get().current_benchmark;
          Hyrise::get().cardinality_statistics += std::to_string(Hyrise::get().current_operator);
          Hyrise::get().cardinality_statistics += ",Aggregate,";
          Hyrise::get().cardinality_statistics += std::to_string(get_table_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(aggregate_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node)) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node)) + "\n";

          // std::cout << "Information for " << *node->lqp_node << std::endl;
          // const auto aggregate_op = std::dynamic_pointer_cast<const AggregateHash>(node);
          // const auto& aggregate_performance_data = aggregate_op->performance_data;
          // const auto& get_table_performance_data = node->left_input()->performance_data;
          // std::cout << "Input size to aggregate: " << get_table_performance_data->output_row_count << std::endl;
          // std::cout << "Result size after aggregate: " << aggregate_performance_data->output_row_count << std::endl;
          // // std::cout << "Actual selectivity of aggregate: " << static_cast<float>(get_table_performance_data->output_row_count) / static_cast<float>(aggregate_performance_data->output_row_count) << std::endl;

          // std::cout << "Input size to aggregate: " << cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node) << std::endl;
          // std::cout << "Result size after aggregate: " << cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node) << std::endl;
        }
      } else if (node->type () == OperatorType::JoinHash) {
        if (node->left_input()->type() == OperatorType::GetTable || node->left_input()->type() == OperatorType::Validate) {
          const auto aggregate_op = std::dynamic_pointer_cast<const JoinHash>(node);
          const auto& aggregate_performance_data = aggregate_op->performance_data;
          const auto& get_table_performance_data = node->left_input()->performance_data;

          Hyrise::get().cardinality_statistics += Hyrise::get().current_benchmark;
          Hyrise::get().cardinality_statistics += std::to_string(Hyrise::get().current_operator);
          Hyrise::get().cardinality_statistics += ",JoinHashLeft,";
          Hyrise::get().cardinality_statistics += std::to_string(get_table_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(aggregate_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node)) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node)) + "\n";

          // std::cout << "Information for " << *node->lqp_node << std::endl;
          // const auto aggregate_op = std::dynamic_pointer_cast<const JoinHash>(node);
          // const auto& aggregate_performance_data = aggregate_op->performance_data;
          // const auto& get_table_performance_data = node->left_input()->performance_data;
          // std::cout << "Left Input size to join: " << get_table_performance_data->output_row_count << std::endl;
          // std::cout << "Result size after join: " << aggregate_performance_data->output_row_count << std::endl;
          // // std::cout << "Actual selectivity of join: " << static_cast<float>(get_table_performance_data->output_row_count) / static_cast<float>(aggregate_performance_data->output_row_count) << std::endl;

          // std::cout << "Left Input size to join: " << cardinality_estimator.estimate_cardinality(node->left_input()->lqp_node) << std::endl;
          // std::cout << "Result size after join: " << cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node) << std::endl;
        }
        if (node->right_input()->type() == OperatorType::GetTable || node->right_input()->type() == OperatorType::Validate) {
          const auto aggregate_op = std::dynamic_pointer_cast<const JoinHash>(node);
          const auto& aggregate_performance_data = aggregate_op->performance_data;
          const auto& get_table_performance_data = node->right_input()->performance_data;

          Hyrise::get().cardinality_statistics += Hyrise::get().current_benchmark;
          Hyrise::get().cardinality_statistics += std::to_string(Hyrise::get().current_operator);
          Hyrise::get().cardinality_statistics += ",JoinHashRight,";
          Hyrise::get().cardinality_statistics += std::to_string(get_table_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(aggregate_performance_data->output_row_count) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(node->right_input()->lqp_node)) + ",";
          Hyrise::get().cardinality_statistics += std::to_string(cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node)) + "\n";

          // std::cout << "Information for " << *node->lqp_node << std::endl;
          // const auto aggregate_op = std::dynamic_pointer_cast<const JoinHash>(node);
          // const auto& aggregate_performance_data = aggregate_op->performance_data;
          // const auto& get_table_performance_data = node->right_input()->performance_data;
          // std::cout << "Right Input size to join: " << get_table_performance_data->output_row_count << std::endl;
          // std::cout << "Result size after join: " << aggregate_performance_data->output_row_count << std::endl;
          // // std::cout << "Actual selectivity of join: " << static_cast<float>(get_table_performance_data->output_row_count) / static_cast<float>(aggregate_performance_data->output_row_count) << std::endl;

          // std::cout << "Right Input size to join: " << cardinality_estimator.estimate_cardinality(node->right_input()->lqp_node) << std::endl;
          // std::cout << "Result size after join: " << cardinality_estimator.estimate_cardinality(aggregate_op->lqp_node) << std::endl;
        }
      }
      
      ++Hyrise::get().current_operator;
      return PQPVisitation::VisitInputs;
    };

    visit_pqp(plan_root, visitor);
  }

  // #####################
  // #####################
  // #####################

  if (pipeline_status == SQLPipelineStatus::Failure) {
    return {pipeline_status, nullptr};
  }
  Assert(pipeline_status == SQLPipelineStatus::Success, "Unexpected pipeline status");

  metrics.emplace_back(std::move(pipeline.metrics()));

  if (expected_result_table) {
    _compare_tables(result_table, expected_result_table, "Using dedicated expected result table");
  } else if (_sqlite_connection) {
    _verify_with_sqlite(pipeline);
  }

  if (_visualize_prefix) {
    _visualize(pipeline);
  }

  return {pipeline_status, result_table};
}

BenchmarkSQLExecutor::~BenchmarkSQLExecutor() {
  if (transaction_context) {
    Assert(transaction_context->phase() == TransactionPhase::Committed ||
               transaction_context->phase() == TransactionPhase::RolledBackByUser ||
               transaction_context->phase() == TransactionPhase::RolledBackAfterConflict,
           "Explicitly created transaction context should have been explicitly committed or rolled back");
  }

  // If the benchmark item does not explicitly manage the transaction life time by calling commit or rollback,
  // we automatically commit the sqlite transaction so that the next item can start a new one:
  if (_sqlite_transaction_open) {
    _sqlite_connection->raw_execute_query("COMMIT TRANSACTION");
  }
}

void BenchmarkSQLExecutor::commit() {
  Assert(transaction_context && !transaction_context->is_auto_commit(),
         "Can only explicitly commit transaction if auto-commit is disabled");
  Assert(transaction_context->phase() == TransactionPhase::Active, "Expected transaction to be active");
  transaction_context->commit();
  if (_sqlite_connection) {
    _sqlite_transaction_open = false;
    _sqlite_connection->raw_execute_query("COMMIT TRANSACTION");
  }
}

void BenchmarkSQLExecutor::rollback() {
  Assert(transaction_context && !transaction_context->is_auto_commit(),
         "Can only explicitly roll back transaction if auto-commit is disabled");
  Assert(transaction_context->phase() == TransactionPhase::Active, "Expected transaction to be active");
  transaction_context->rollback(RollbackReason::User);
  if (_sqlite_connection) {
    _sqlite_transaction_open = false;
    _sqlite_connection->raw_execute_query("ROLLBACK TRANSACTION");
  }
}

void BenchmarkSQLExecutor::_verify_with_sqlite(SQLPipeline& pipeline) {
  Assert(pipeline.statement_count() == 1, "Expecting single statement for SQLite verification");

  const auto sqlite_result = _sqlite_connection->execute_query(pipeline.get_sql());
  const auto [pipeline_status, result_table] = pipeline.get_result_table();
  Assert(pipeline_status == SQLPipelineStatus::Success, "Non-successful pipeline should have been caught earlier");

  // Modifications (INSERT, UPDATE, DELETE) do not return a table. We do not know what changed - we do not even know
  // which table has been modified. Extracting that info from the plan and verifying the entire table would take way
  // too long. As such, we rely on any errors here to be found when the potentially corrupted data is SELECTed the
  // next time.
  if (!result_table) {
    return;
  }

  _compare_tables(result_table, sqlite_result, "Using SQLite's result table as expected result table");
}

void BenchmarkSQLExecutor::_compare_tables(const std::shared_ptr<const Table>& actual_result_table,
                                           const std::shared_ptr<const Table>& expected_result_table,
                                           const std::optional<const std::string>& description) {
  Timer timer;

  if (actual_result_table->row_count() > 0) {
    Assert(expected_result_table, "Verifying SQL statement on sqlite failed: No SQLite result");
    if (expected_result_table->row_count() == 0) {
      any_verification_failed = true;
      if (description) {
        std::cout << "- " + *description << "\n";
      }
      std::cout << "- Verification failed: Hyrise's actual result is not empty, but the expected result is ("
                << timer.lap_formatted() << ")"
                << "\n";
    } else if (const auto table_difference_message = check_table_equal(
                   actual_result_table, expected_result_table, OrderSensitivity::No, TypeCmpMode::Lenient,
                   FloatComparisonMode::RelativeDifference, IgnoreNullable::Yes)) {
      any_verification_failed = true;
      if (description) {
        std::cout << *description << "\n";
      }
      std::cout << "- Verification failed (" << timer.lap_formatted() << ")"
                << "\n"
                << *table_difference_message << "\n";
    }
  } else {
    if (expected_result_table && expected_result_table->row_count() > 0) {
      any_verification_failed = true;
      if (description) {
        std::cout << *description << "\n";
      }
      std::cout << "- Verification failed: Expected result table is not empty, but Hyrise's actual result is ("
                << timer.lap_formatted() << ")"
                << "\n";
    }
  }
}

void BenchmarkSQLExecutor::_visualize(SQLPipeline& pipeline) {
  GraphvizConfig graphviz_config;
  graphviz_config.format = "svg";

  const auto& lqps = pipeline.get_optimized_logical_plans();
  const auto& pqps = pipeline.get_physical_plans();

  auto prefix = *_visualize_prefix;

  if (_num_visualized_plans == 1) {
    // We have already visualized a prior SQL pipeline in this benchmark item - rename the existing file
    std::filesystem::rename(prefix + "-LQP.svg", prefix + "-0-LQP.svg");
    std::filesystem::rename(prefix + "-PQP.svg", prefix + "-0-PQP.svg");
  }
  if (_num_visualized_plans > 0) {
    prefix += "-" + std::to_string(_num_visualized_plans);
  }

  LQPVisualizer{graphviz_config, {}, {}, {}}.visualize(lqps, prefix + "-LQP.svg");
  PQPVisualizer{graphviz_config, {}, {}, {}}.visualize(pqps, prefix + "-PQP.svg");

  ++_num_visualized_plans;
}

}  // namespace hyrise
