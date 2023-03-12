#!/usr/bin/env bash

histogram="MaxDiffAreaHistogram"

curr_date=$(date +"%Y-%m-%d-%T")
thread_count=$1

repitions=1
scale_factor=100

cd ../cmake-build-release/
ninja 
cd ../scripts

setup_env() {
  histogram=$1
  benchmark_name=$2

  build_times_filename="build_times_multithreaded_$curr_date.csv"

  echo "Measuring $benchmark_name with $histogram ..."
  echo "Thread Count: $thread_count"
  export HISTOGRAM="$histogram"
  export BUILD_TIME="$build_times_filename"
  export THREADCOUNT="$thread_count"
}

get_tpcds_data() {
  benchmark_name="tpcds"
  setup_env $1 $benchmark_name

  ./cmake-build-release/hyriseBenchmarkTPCDS -r $repitions -s $scale_factor --scheduler
}


get_tpcds_data $histogram
