#!/usr/bin/env bash

###########################################################################
#
# Builds TPCH with all available Histograms and fetches the data.
#
# Dependencies of this script:
#     1. Execute in scripts folder (I was too lazy to remove that dependency).
#     2. The build directory needs to called "cmake-build-debug".
#     3. The build of Hyrise needs to be done and finished.
#     4. The jupyter data needs to be stored in jupyter/data/.
#
# Resulting data will be of name <executed-binary>-<histogram>-cardinalities-<timestamp>.csv
#
###########################################################################

histograms=("EquiHeightHistogram" "EquiWidthHistogram" "EqualDistinctCountHistogram" "MaxDiffFrequencyHistogram")

if [[ $# -lt 2 ]]; then
  echo "Usage: ./get_histogram_data.sh benchmark_binary args (dont forget to put brackets around the args"
  exit 1
fi

build_times_filename="build_times.csv"

get_histogram_data () {
  curr_date=$(date +"%Y-%m-%d-%T")
  echo $curr_date

  args="${*:3}"
  histogram=$2
  binary=$1

  cardinalities_filename="$binary-$histogram-cardinalities-$curr_date.csv"

  echo "Measuring $histogram ($binary $args) ..."
  export HISTOGRAM="$histogram"
  export CARDINALITIES="$cardinalities_filename"
  export BUILD_TIME="$build_times_filename"

  cd ../cmake-build-debug/
  #ninja

  ./$binary $args

  mv $cardinalities_filename ../jupyter/data
}

for hist in ${histograms[@]}; do
  get_histogram_data $1 $hist $2
done

mv $build_times_filename ../jupyter/data

