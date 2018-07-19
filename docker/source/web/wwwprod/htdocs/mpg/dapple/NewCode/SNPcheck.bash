#!/bin/bash

awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr1.ws > SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr2.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr3.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr4.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr5.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr6.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr7.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr8.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr9.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr10.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr11.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr12.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr13.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr14.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr15.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr16.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr17.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr18.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr19.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr20.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr21.ws >> SNPcheck
awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/chr22.ws >> SNPcheck

