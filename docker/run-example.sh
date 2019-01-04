#!/usr/bin/env bash
set -euo pipefail

#
# Example dapple command
#
#   This script demonstrates how to run the 
# dapple command in a docker container
#


#
# 'docker_img', default: dapple:0.18.3
#
: ${docker_img=genepattern/dapple:v0.19-pre.4}

#
# 'dapple_dir', default: ../
#
: ${dapple_dir=$(dirname `pwd`)}

#
# 'working_dir', default: ./tmp/exampleOutput
#
: ${working_dir=`pwd`/tmp/exampleOutput}
mkdir -p ${working_dir}
rm -f ${working_dir}/*


#
# 'dapple_script', default: /code/BuildNetwork_slowTMP.py
#
# DOCKER
#: ${dapple_script=/code/BuildNetwork_slowTMP.py}
# DEBUG
: ${dapple_script=${dapple_dir}/docker/source/code/BuildNetwork_slowTMP.py}


#
# 'input_file', default: ../exampleInput2
#
#: ${input_file=${dapple_dir}/gpunit/input/SNP_input}
: ${input_file=${dapple_dir}/gpunit/input/GeneNames_input}

echo "Running dapple command ..."
echo "  working_dir: ${working_dir}"

set -x

docker run \
  -w "${working_dir}" \
  -v "${dapple_dir}":/"${dapple_dir}" \
  ${docker_img} \
python "${dapple_script}" \
  "${input_file}" \
  keyword=SNPpermuteHG18 \
  permute=50 \
  plot=false \
  genome=18 \
  nearestgene=false \
  seed=123

set +x

exit_code=$?
echo "Completed, see results in ${working_dir}"

exit $exit_code
