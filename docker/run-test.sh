#!/usr/bin/env bash
set -euo pipefail

# Example command line
# docker run \
#   -w `pwd`/tmp/jobResults/1 \
#   -v "$(dirname `pwd`)":"$(dirname `pwd`)" \
#   -v `pwd`/source/web:/web \
#   -v `pwd`/source/fg:/fg \
#   -v `pwd`/source/home:/home \
#   dapple:0.18.3 \
#     python \
#     /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py \
#     "$(dirname `pwd`)"/exampleInput \
#     keyword=example \
#     permute=0 \
#     plot=true \
#     genome=19 \
#     nearestgene=false \
# >`pwd`/tmp/jobResults/1/stdout.txt 2>`pwd`/tmp/jobResults/1/stderr.txt


# Paramaterized command line

# set the docker image
: ${docker_img=dapple:0.18.3}

# set the root directory 
: ${dapple_dir=$(dirname `pwd`)}

# set the working directory, default=tmp/jobResults/1
#   To match the way GenePattern Server runs jobs
#
: ${job_id=1}
: ${working_dir=`pwd`/tmp/jobResults/${job_id}}

mkdir -p ${working_dir}

# set the source dir
: ${source_dir=${dapple_dir}/docker/source}

# set the input file
: ${input_file=${dapple_dir}/exampleInput}

docker run \
  -w "${working_dir}" \
  -v "${dapple_dir}":"${dapple_dir}" \
  -v "${source_dir}/web/wwwprod/htdocs/mpg/dapple:/web/wwwprod/htdocs/mpg/dapple" \
  -v "${source_dir}/fg:/fg" \
  -v "${source_dir}/home:/home" \
  ${docker_img} \
    python \
    /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py \
    "${input_file}" \
    keyword=example \
    permute=0 \
    plot=true \
    genome=19 \
    nearestgene=false \
>${working_dir}/stdout.txt 2>${working_dir}/stderr.txt
