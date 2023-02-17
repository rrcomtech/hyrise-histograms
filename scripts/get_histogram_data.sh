#!/usr/bin/env bash

histograms=("EquiHeightHistogram" "EquiWidthHistogram" "EqualDistinctCountHistogram" "MaxDiffFrequencyHistogram" "GDYHistogram")

build_times_filename="build_times.csv"
curr_date=$(date +"%Y-%m-%d-%T")

repitions=1
scale_factor=10

results_folder="$(pwd)/../cmake-build-debug/results"
mkdir $results_folder

cd ../cmake-build-debug/
ninja 
cd ../scripts

setup_env() {
  histogram=$1
  benchmark_name=$2

  cardinalities_filename="$benchmark_name-$histogram-cardinalities-$curr_date.csv"
  extended_build_times_filename="$benchmark_name-$build_times_filename"

  echo "Measuring $benchmark_name with $histogram ..."
  export HISTOGRAM="$histogram"
  export CARDINALITIES="$cardinalities_filename"
  export BUILD_TIME="$extended_build_times_filename"

  mkdir "$results_folder/$benchmark_name/"
}

get_tpch_data() {
  benchmark_name="tpch"
  setup_env $1 $benchmark_name

  cd ../cmake-build-debug/

  ./hyriseBenchmarkTPCH -r $repitions -s $scale_factor

  mv $cardinalities_filename "$results_folder/$benchmark_name/"
  mv $extended_build_times_filename "$results_folder/$benchmark_name/"
}

get_job_data() {
  benchmark_name="job"
  setup_env $1 $benchmark_name

  cd ../

  ./cmake-build-debug/hyriseBenchmarkJoinOrder

  mv $cardinalities_filename "$results_folder/$benchmark_name/"
  mv $extended_build_times_filename "$results_folder/$benchmark_name/"

  # Just to be in an expected state for the next benchmark.
  cd scripts
}

get_tpcds_data() {
  benchmark_name="tpcds"
  setup_env $1 $benchmark_name

  cd ../cmake-build-debug/

  ./hyriseBenchmarkTPCDS -r $repitions -s $scale_factor

  mv $cardinalities_filename "$results_folder/$benchmark_name/"
  mv $extended_build_times_filename "$results_folder/$benchmark_name/"
}

get_tpcc_data() {
  benchmark_name="tpcc"
  setup_env $1 $benchmark_name

  cd ../cmake-build-debug/

  ./hyriseBenchmarkTPCC -r $repitions -s $scale_factor

  mv $cardinalities_filename "$results_folder/$benchmark_name/"
  mv $extended_build_times_filename "$results_folder/$benchmark_name/"
}

for hist in ${histograms[@]}; do
  get_tpch_data $hist
  get_job_data $hist
  get_tpcds_data $hist
  get_tpcc_data $hist
done