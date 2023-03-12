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

thread_count=1
while [ $thread_count -lt 65 ]
do
  get_tpcds_data $histogram $thread_count
  let thread_count=thread_count+1
done
