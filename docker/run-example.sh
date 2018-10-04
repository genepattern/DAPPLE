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
: ${docker_img=dapple:0.19}

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
# 'source_dir', default: ./source
#
#: ${source_dir=${dapple_dir}/docker/source}

#
# 'dapple_script', default: /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py
#
: ${dapple_script=/code/BuildNetwork_slowTMP.py}


#
# 'input_file', default: ../exampleInput
#
#: ${input_file=${dapple_dir}/exampleInput}
: ${input_file=${dapple_dir}/exampleInput2}
#: ${input_file=${dapple_dir}/exampleInput_SNP_RA}
#: ${input_file=${dapple_dir}/exampleInput_SNP_CD}


echo "Running dapple command ..."
echo "  working_dir: ${working_dir}"

set -x

docker run \
  -w "${working_dir}" \
  -v "${dapple_dir}":/"${dapple_dir}" \
  ${docker_img} \
python "${dapple_script}" \
  "${input_file}" \
  keyword=p1000NearestFalseHG18try2 \
  permute=1000 \
  plot=true \
  genome=18 \
  nearestgene=false 

#  -v "${source_dir}/code":/code \

  #/bin/bash -c "awk 'NR==FNR{s[\$1]=1;next}{if (s[\$1]==1) print }' Seeds <(gzip -dc /data/wingspan/hg18/*.ws.gz |) > SNPcheck"
#head <(/bin/bash -c "gzip -dc /data/wingspan/hg18/*.ws.gz |")
#awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds `/bin/bash -c "ls /data/wingspan/hg18/*.ws.gz | xargs gzip -dc"`
#awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds <(/bin/bash -c "gzip -dc /data/wingspan/hg18/*.ws.gz")

#ls /data/wingspan/hg18/chr1.ws.gz

#  /bin/bash -c "ls /data/wingspan/hg18/*"

#  -v "${source_dir}/data":/data \

#  -v "${source_dir}/web/wwwprod/htdocs/mpg/dapple/NewCode":/web/wwwprod/htdocs/mpg/dapple/NewCode \
#  -v "${source_dir}/fg/wgas2/rossin":/fg/wgas2/rossin \
#  -v "${source_dir}/home/unix/rossin/DALY/PPI/NewCode":/home/unix/rossin/DALY/PPI/NewCode \

set +x

# >${working_dir}/stdout.txt 2>${working_dir}/stderr.txt
exit_code=$?
echo "Completed, see results in ${working_dir}"

exit $exit_code
