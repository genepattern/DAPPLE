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
: ${docker_img=dapple:0.18.3}

#
# 'dapple_dir', default: ../
#
: ${dapple_dir=$(dirname `pwd`)}

#
# 'working_dir', default: ./tmp/exampleOutput
#
: ${working_dir=`pwd`/tmp/exampleOutput}
mkdir -p ${working_dir}

#
# 'source_dir', default: ./source
#
: ${source_dir=${dapple_dir}/docker/source}

#
# 'dapple_script', default: /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py
#
: ${dapple_script=/web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py}


#
# 'input_file', default: ../exampleInput
#
: ${input_file=${dapple_dir}/exampleInput}

echo "Running dapple command ..."
echo "  working_dir: ${working_dir}"

docker run \
  -w "${working_dir}" \
  -v "${dapple_dir}":"${dapple_dir}" \
  -v "${source_dir}/web/wwwprod/htdocs/mpg/dapple/NewCode":/web/wwwprod/htdocs/mpg/dapple/NewCode \
  -v "${source_dir}/fg/wgas2/rossin":/fg/wgas2/rossin \
  -v "${source_dir}/home/unix/rossin/DALY/PPI/NewCode":/home/unix/rossin/DALY/PPI/NewCode \
  ${docker_img} \
python "${dapple_script}" \
  "${input_file}" \
  keyword=example \
  permute=0 \
  plot=true \
  genome=19 \
  nearestgene=false 

# >${working_dir}/stdout.txt 2>${working_dir}/stderr.txt
exit_code=$?
echo "Completed, see results in ${working_dir}"

exit $exit_code
