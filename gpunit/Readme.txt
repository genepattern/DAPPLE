Three test:

Set seed to 123
Set iterations to 50

Fast [15s]:
Am I alive. Only uses gene names, no permutation
Check if DAPPLE can a) load a network and b) generate simple output

GeneNamesPermute [200s]:
Uses gene names, permutes to calculate p values. Checks for a), b) and c) that permutation works

SNP_permute [550s]
Uses SNPs as input, permutes to calculate p values. Checks a), b), c), and d) SNPs can be mapped to regions/genes

TODO:
Check all three genome versions [18, 19, 1kg]

