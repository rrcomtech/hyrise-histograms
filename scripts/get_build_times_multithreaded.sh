#!/usr/bin/env bash

histogram="MaxDiffAreaHistogram"

curr_date=$(date +"%Y-%m-%d-%T")

cd ./cmake-build-release/
ninja 
cd ..

setup_env() {
  histogram=$1
  benchmark_name=$2
  thread_count=$3

  build_times_filename="build_times_multithreaded_$curr_date.csv"

  echo "Measuring $benchmark_name with $histogram ..."
  echo "Thread Count: $thread_count"
  export HISTOGRAM="$histogram"
  export BUILD_TIME="$build_times_filename"
  export THREADCOUNT="$thread_count"
}

get_tpcds_data() {
  benchmark_name="tpcds"
  setup_env $1 $benchmark_name $2

  ./cmake-build-release/hyriseBenchmarkTPCDS -r 1 -s 100 --scheduler
}

get_tpcds_data $histogram 1
get_tpcds_data $histogram 2
get_tpcds_data $histogram 3
get_tpcds_data $histogram 4
get_tpcds_data $histogram 6

thread_count=8
while [ $thread_count -lt 65 ]
do
  get_tpcds_data $histogram $thread_count
  let thread_count=thread_count+4
done
