#!/bin/bash

URL="https://personal.broadinstitute.org/hhorn/dappleData/"

mkdir -p /data/bim
mkdir -p /data/bim/hg18
mkdir -p /data/bim/hg19
mkdir -p /data/bim/1KG

mkdir -p /data/wingspan
mkdir -p /data/wingspan/hg18
mkdir -p /data/wingspan/hg19
mkdir -p /data/wingspan/1KG

for chr in {1..22}; do
    wget ${URL}bim/1KG/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1_chr${chr}.eur.bfile.bim.gz -q -O /data/bim/1KG/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1_chr${chr}.eur.bfile.bim.gz
    wget ${URL}bim/hg18/chr${chr}.BIM.gz -q -O /data/bim/hg18/chr${chr}.BIM.gz
    wget ${URL}bim/hg19/chr${chr}.BIM.gz -q -O /data/bim/hg19/chr${chr}.BIM.gz

    wget ${URL}wingspan/hg18/chr${chr}.ws.gz -q -O /data/wingspan/hg18/chr${chr}.ws.gz
    wget ${URL}wingspan/hg19/chr${chr}.ws.gz -q -O /data/wingspan/hg19/chr${chr}.ws.gz
    wget ${URL}wingspan/1KG/chr${chr}.ws.gz -q -O /data/wingspan/1KG/chr${chr}.ws.gz
done
