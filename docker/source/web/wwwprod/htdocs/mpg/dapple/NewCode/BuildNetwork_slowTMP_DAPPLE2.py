import os
import sys
import re
import random
import copy
import json
import numpy
from Functions_slowTMP_DAPPLE2 import *

##input can be SNPs, genes or regions

inputlist = open(sys.argv[1],'r').readlines()
inputlist = list(set(inputlist))
genestospecify = "X"
CIcutoff = 2
permutlength=1000
plot = False
plotP = False
manualColoring = False
iterate = False
compare = False
regulatory = [50000,50000]
dir = "./"
nearestgene = False
seedRandom = "NA"
collapseCI = False
genomebuild="19"
def str2bool(v):
	return v.lower() in ("true")
if len(sys.argv)>2:
	for argv in sys.argv[2:]:
		if argv.split("=")[0] == "specify":
			genestospecify = open(argv.split("=")[1],'r')
		elif argv.split("=")[0] == "permute":
			permutlength = int(argv.split("=")[1])
		elif argv.split("=")[0] == "CIcutoff":
			CIcutoff = int(argv.split("=")[1])
		elif argv.split("=")[0] == "plot":
			plot = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "iterate":
			iterate = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "manual":
			manualColoring = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "keyword":
			sys.argv[1] = argv.split("=")[1]
		elif argv.split("=")[0] == "compare":
			compare = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "regulatory":
			regulatory = [1000*int(argv.split("=")[1].split(",")[0]),1000*int(argv.split("=")[1].split(",")[1])]
		elif argv.split("=")[0] == "dir":
			dir = argv.split("=")[1]
		elif argv.split("=")[0] == "nearestgene":
			nearestgene = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "seed":
			seedRandom = int(argv.split("=")[1])
		elif argv.split("=")[0] == "plotP":
			plotP = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "collapseCI":
			collapseCI = str2bool(argv.split("=")[1])
		elif argv.split("=")[0] == "genome":
			genomebuild = argv.split("=")[1]


##Change to the specified working directory
os.chdir(dir)

##Create gene name conversion dictionary
genenames = open("/home/radon00/rossin/DALY/PPI/NewCode/IWtoHugo",'r').readlines()
geneNameConvert={}
for name in genenames:
	name = name.strip("\n").split("\t")
	geneNameConvert[int(name[0])] = name[1]
	geneNameConvert[name[1]]=int(name[0])

##load PPI data and create dictionaries
ppifile = "/fg/wgas2/rossin/KasperData/InWeb3_HC_Red_sorted"

##create a list of proteins in the high-confidence set
hcProteins = []
ppi = open(ppifile,'r')
for line in ppi:
	line = line.strip("\n").split("\t")
	hcProteins.append(int(line[0]))
	hcProteins.append(int(line[1]))
hcProteins = list(set(hcProteins))

##set gene input as false
geneInput=False
snps = []
seeds = []

##Open Missing Genes file
outfile=open(sys.argv[1]+"_MissingGenes",'w')
outfile.write("SNP\tNumGenes\tNumGenes_InDatabase\tPercent\tGenesInDatabase\tMissingGenes\n")
outfile.close()

##if input is a list of SNPs and nearestgene = true
matched = []
if nearestgene and anySnp(inputlist):
	print("Input has snps and you have chosen to find the nearest gene")
	snpList = []
	for line in inputlist:
		if isSnp(line):
			snpList.append(line.strip("\n").strip("\r"))
			matched.append(inputlist.index(line))
	x=NearestGene(snpList,build=genomebuild)
	seedsNG = x[0]
	for i in range(len(seedsNG)):
		seedsNG[i]=geneNameConvert[seedsNG[i]]
	seeds = seeds+seedsNG
	snps = snps+x[1]

##if input is a list of SNPs
if anySnp(inputlist) and not nearestgene:
	print("Input has snps. Converting snps to genes...")
	snpList=[]
	for line in inputlist:
		if isSnp(line):
			snpList.append(line.strip("\n").strip("\r"))
			matched.append(inputlist.index(line))
	x = SNPtoGene(snpList,regulatory[0],regulatory[1],build=genomebuild)
	seeds = seeds+x[0]
	for i in range(len(seeds)):
		seeds[i]=geneNameConvert[seeds[i]]
	snps = snps+x[1]

##if input is a list of genes
if anyGene(inputlist) or anyGeneRegion(inputlist):
	print("Input has genes. Extracting genes...")
	inputlist = list(set(inputlist))
	geneInput = True
	count=0
	missingSeeds=[]
	missingSNPs=[]
	for line in inputlist:
		if not isGene(line) and not isGeneRegion(line):
			continue
		matched.append(inputlist.index(line))
		if geneNameConvert.has_key(line.strip("\n").split()[0]):
			if geneNameConvert[line.strip("\n").split()[0]] in hcProteins:
				seeds.append(geneNameConvert[line.strip("\n").split()[0]])
				if len(inputlist[0].split())==2:
					snps.append(line.strip("\n").split()[1])
				else:
					snps.append("G"+str(count))
					count+=1
			else:
				missingSeeds.append(geneNameConvert[line.strip("\n").split()[0]])
				if len(inputlist[0].split())==2:
					missingSNPs.append(line.strip("\n").split()[1])
				else:
					missingSNPs.append("gene"+str(count))
					count+=1

##if input is a list of regions

if anyRegion(inputlist):
	print("Input has regions. Converting regions to genes...")
	for line in inputlist:
		if isRegion(line):
			matched.append(inputlist.index(line))
	x = RegionToGene(inputlist,regulatory[0],regulatory[1],build=genomebuild)
	seedsRegion = x[0]
	for i in range(len(seedsRegion)):
		seedsRegion[i]=geneNameConvert[seedsRegion[i]]
	seeds=seeds+seedsRegion
	snps = snps+x[1]

unmatched = set(range(0,len(inputlist)))-set(matched)
errorInputs = []
for i in unmatched:
	errorInputs.append(inputlist[i])

if len(unmatched)>0:
	print("The following inputs are invalid:\n"+"\n".join(errorInputs))

##Exit if too few genes match PPI data
if float(len(seeds))/float(len(inputlist)) < 0.1:
	print("Less than 10% of seed genes are in the PPI database. That is too few. Check to make sure IDs are Hugo IDs - for example, \"ATXN1\". If your input was a list of SNPs, make sure they are in HapMap release 22.")
	sys.exit()
	
##merge overlapping regions or genes

merged = mergeRegions(seeds,snps)
seeds=merged[0]
snps=merged[1]


##create gene list

if genestospecify !="X":
	tokeep = []
	for line in genestospecify:
		line=line.rstrip()
		if geneNameConvert.has_key(line):
			tokeep.append(geneNameConvert[line.strip("\n")])
else:
	tokeep = []

tokeepOriginal = tokeep[:]

SNPtoRemove = []
for i in range(len(seeds)):
	if seeds[i] in tokeep:
		SNPtoRemove.append(snps[i])
indicestoremove = []
for i in range(len(snps)):
	if snps[i] in SNPtoRemove:
		if seeds[i] in tokeep:
			continue
		else:
			indicestoremove.append(i)

for index in indicestoremove:
	seeds[index] = 'X'
	snps[index] = 'X'
while 'X' in seeds:
	seeds.remove('X')
while 'X' in snps:
	snps.remove('X')


genesnp = {}
snpGeneCount = {}
for i in range(len(seeds)):	
	genesnp[seeds[i]] = snps[i]
	if not snpGeneCount.has_key(snps[i]):
		snpGeneCount[snps[i]]=0
	snpGeneCount[snps[i]]+=1

genesnpCopy = copy.deepcopy(genesnp)

originalseeds=seeds[:]
originalsnps=snps[:]

##write out missing genes file if geneInput
if (geneInput):
	missingGeneSnp = {}
	missingGeneCount = {}
	for i in range(len(missingSeeds)):
		missingGeneSnp[missingSeeds[i]] = missingSNPs[i]
		if not missingGeneCount.has_key(missingSNPs[i]):
			missingGeneCount[missingSNPs[i]]=0
		missingGeneCount[missingSNPs[i]]+=1
	for snp in snps:
		if not missingGeneCount.has_key(snp):
			missingGeneCount[snp]=0
	for snp in missingSNPs:
		if not snpGeneCount.has_key(snp):
			snpGeneCount[snp]=0
	outfile = open(sys.argv[1]+"_MissingGenes",'a')
	allSnps = list(set(snps+missingSNPs))
	for snp in allSnps:
		outfile.write(snp+"\t"+str(snpGeneCount[snp]+missingGeneCount[snp])+"\t"+str(snpGeneCount[snp])+"\t")
		outfile.write(str(round(float(snpGeneCount[snp])/int(snpGeneCount[snp]+missingGeneCount[snp]),2))+"\t")
		for gene in genesnp.keys():
			if genesnp[gene]==snp:
				outfile.write(geneNameConvert[gene]+",")
		if not snp in genesnp.values():
			outfile.write("-")
		outfile.write("\t")
		for gene in missingGeneSnp.keys():
			if missingGeneSnp[gene]==snp:
				outfile.write(geneNameConvert[gene]+",")
		if not snp in missingGeneSnp.keys():
			outfile.write("-")
		outfile.write("\n")

outfile.close()

##Organize missing gene file
#os.system("sort "+sys.argv[1]+"_MissingGenes | uniq > "+sys.argv[1]+"_MissingGenesTmp")
#os.system("mv "+sys.argv[1]+"_MissingGenesTmp "+sys.argv[1]+"_MissingGenes")

##print seed-locus file
seedLocusFile = open(sys.argv[1]+"_seedLocusMapping",'w')
for seed in seeds:
	seedLocusFile.write(geneNameConvert[seed]+"\t"+genesnp[seed]+"\n")
seedLocusFile.close()

##run disease network, set up dictionaries
print("Building disease network...")

diseaseNet = buildNet(ppifile,"/home/radon00/rossin/DALY/PPI/NewCode/ppiannot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])


seedSeed = diseaseNet["seedSeed"]
seedCI = diseaseNet["seedCI"]
CIseed = diseaseNet["CIseed"]
CIdegreesDict = diseaseNet["CIdegreesDict"]
CIdegreesMean = diseaseNet["CIdegreesMean"]
geneCountDirect = diseaseNet["geneCountDirect"]
geneCountIndirect = diseaseNet["geneCountIndirect"]
locusCountDirect = diseaseNet["locusCountDirect"]
locusCountIndirect = diseaseNet["locusCountIndirect"]
gene2locusEdgeCountDirect = diseaseNet["gene2locusEdgeCountDirect"]
gene2locusEdgeCountIndirect = diseaseNet["gene2locusEdgeCountIndirect"]
locus2locusEdgeCountDirect = diseaseNet["locus2locusEdgeCountDirect"]
locus2locusEdgeCountIndirect = diseaseNet["locus2locusEdgeCountIndirect"]
gene2locusDegreeDictDirect = diseaseNet["gene2locusDegreeDictDirect"]
gene2locusDegreesMeanDirect = diseaseNet["gene2locusDegreesMeanDirect"]
gene2locusDegreeDictIndirect = diseaseNet["gene2locusDegreeDictIndirect"]
gene2locusDegreesMeanIndirect = diseaseNet["gene2locusDegreesMeanIndirect"]
gene2locusWeightedDegreeDictIndirect = diseaseNet["gene2locusWeightedDegreeDictIndirect"]
gene2locusWeightedDegreesMeanIndirect = diseaseNet["gene2locusWeightedDegreesMeanIndirect"]
locus2locusDegreeDictDirect = diseaseNet["locus2locusDegreeDictDirect"]
locus2locusDegreesMeanDirect = diseaseNet["locus2locusDegreesMeanDirect"]
locus2locusDegreeDictIndirect = diseaseNet["locus2locusDegreeDictIndirect"]
locus2locusDegreesMeanIndirect = diseaseNet["locus2locusDegreesMeanIndirect"]
locus2locusWeightedDegreeDictIndirect = diseaseNet["locus2locusWeightedDegreeDictIndirect"]
locus2locusWeightedDegreesMeanIndirect = diseaseNet["locus2locusWeightedDegreesMeanIndirect"]


##Write out indirect network files
outfile = open(sys.argv[1]+"_indirectNetwork",'w')
icheck = False
for seed in seedCI.keys():
	if len(seedCI[seed])>0:
		icheck=True
		for CI in seedCI[seed]:
			outfile.write(geneNameConvert[seed]+" -- "+geneNameConvert[CI]+"\n")

outfile.close()

##Write out direct network files
outfile = open(sys.argv[1]+"_directNetwork",'w')
tally = []
dcheck=False
for seed in seedSeed.keys():
	if len(seedSeed[seed]) > 0:
		dcheck=True
		for otherSeed in seedSeed[seed]:
			x = [seed,otherSeed]
			x.sort()
			if not x in tally:
				outfile.write(geneNameConvert[seed]+" -- "+geneNameConvert[otherSeed]+"\n")
				tally.append(x)

outfile.close()

if not icheck and not dcheck:
	plot=False

if (plot and plotP) or (plot and permutlength==0):
	nodemap = {}
	nodes = []
	links = []
	links2 = []
	ct = 0
	plotfile = open('data.json','w')
	plotfile2 = open('dataDirect.json','w')
	directSeedList = []
	if plot:
		print("Plotting network...")
		directOut = open("./tmpDirect",'w')
		indirectOut = open("./tmpIndirect",'w')
		for seed in seedSeed:
			source = geneNameConvert[seed]
			if len(seedSeed[seed])>0:
				directSeedList.append(source)
			for otherSeed in seedSeed[seed]:
				pair = [seed,otherSeed]
				pair.sort()
				directOut.write(geneNameConvert[pair[0]]+"\t"+geneNameConvert[pair[1]]+"\n")
				target = geneNameConvert[otherSeed]
				if not nodemap.has_key(source):
					nodemap[source] = ct
					ct += 1
				if not nodemap.has_key(target):
					nodemap[target] = ct
					ct += 1
				d = {
					"source": nodemap[source],
					"target": nodemap[target],
					"value": 2
					}
				links.append(d)
				links2.append(d)
		directOut.close()
		os.system("sort -k 1,1 tmpDirect | uniq > tmp")
		os.system("mv tmp tmpDirect")
		for seed in seedCI:
			source = geneNameConvert[seed]
			for ci in seedCI[seed]:
				target = geneNameConvert[ci]
				pair = [seed,ci]
				pair.sort()
				indirectOut.write(geneNameConvert[pair[0]]+"\t"+geneNameConvert[pair[1]]+"\n")
				if not nodemap.has_key(source):
					nodemap[source] = ct
					ct += 1
				if not nodemap.has_key(target):
					nodemap[target] = ct
					ct += 1					
				d = {
					"source": nodemap[source],
					"target": nodemap[target],
					"value": 1
					}
				links.append(d)
		indirectOut.close()
		os.system("sort -k 1,1 tmpIndirect | uniq > tmp")
		os.system("mv tmp tmpIndirect")
		outfile = open("tmpSnpGene",'w')
		colormap = {}
		col = 1
		snps = list(set(genesnp.values()))
		for seed in genesnp.keys():
			outfile.write(genesnp[seed]+"\t"+geneNameConvert[seed]+"\n")
			colormap[geneNameConvert[seed]] = int(snps.index(genesnp[seed]))
		outfile.close()
		if not collapseCI:
			collapse = "FALSE"
		else:
			collapse = "TRUE"
		if not manualColoring:
			if not plotP:
				os.system("sed s/manualColoring/FALSE/g /web/wwwprod/htdocs/mpg/dapple/NewCode/plot.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")
			else:
				os.system("sed s/manualColoring/FALSE/g /web/wwwprod/htdocs/mpg/dapple/NewCode/plotP.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")			
		else:
			os.system("sed \'s/manualColoring/TRUE/g\' /web/wwwprod/htdocs/mpg/dapple/NewCode/plot.R | sed \'s/coloringScheme/"+manualColoring+"/g\' | sed \'s/KEYWORDHERE/"+sys.argv[1]+"/g\' > tmpPlot.R")
		os.system("R --vanilla < tmpPlot.R")

##get coordinates
	coordFile = open("PlotCoordinates",'r')
	coords = {}
	for line in coordFile:
		line = line.strip("\n").split()
		coords[line[0]] = [float(line[1]),float(line[2])]

##create data for SVG plot
	if len(seeds)<400:
		nodes2=[]
		for k in nodemap.keys():
			color = colormap.get(k, 0)
			d = {
				"index": nodemap[k],
				"name": k,
				"group": color,
				"x": coords[k][0],
				"y": coords[k][1]
				}
			nodes.append(d)
			if k in directSeedList:
				nodes2.append(d)

		nodes = sorted(nodes, key=lambda x: x['index'])
		nodes2 = sorted(nodes2, key=lambda x: x['index'])
		out2 = {
			"nodes": nodes2,
			"links": links2
			}
		json.dump(out2, plotfile2)
		out = {
			"nodes": nodes,
			"links": links
			}
		json.dump(out, plotfile)
		plotfile.close()
		plotfile2.close()


##exit if permutlength=0
if permutlength==0:
	sys.exit()

##permute
geneCountDirectPermut = ()
geneCountIndirectPermut = ()
locusCountDirectPermut = ()
locusCountIndirectPermut = ()
gene2locusEdgeCountDirectPermut = ()
gene2locusEdgeCountIndirectPermut = ()
locus2locusEdgeCountDirectPermut = ()
locus2locusEdgeCountIndirectPermut = ()
gene2locusDegreesMeanDirectPermut = ()
gene2locusDegreesMeanIndirectPermut = ()
gene2locusWeightedDegreesMeanIndirectPermut = ()
locus2locusDegreesMeanDirectPermut = ()
locus2locusDegreesMeanIndirectPermut = ()
locus2locusWeightedDegreesMeanIndirectPermut = ()
CIdegreesMeanPermut = ()


gene2locusDegreeDictDirectPermut = {}
gene2locusDegreeDictIndirectPermut = {}
gene2locusWeightedDegreeDictIndirectPermut = {}
locus2locusDegreeDictDirectPermut = {}
locus2locusDegreeDictIndirectPermut = {}
locus2locusWeightedDegreeDictIndirectPermut = {}
CIdegreesDictPermut = {}

for seed in seeds:
	gene2locusDegreeDictDirectPermut[seed] = 0
	gene2locusDegreeDictIndirectPermut[seed] = 0
	gene2locusWeightedDegreeDictIndirectPermut[seed] = 0
	locus2locusDegreeDictDirectPermut[seed] = 0
	locus2locusDegreeDictIndirectPermut[seed] = 0
	locus2locusWeightedDegreeDictIndirectPermut[seed] = 0
	CIdegreesDictPermut[seed] = 0

ppifiles = os.listdir("/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/")
for i in range(len(ppifiles)):
	ppifiles[i] = re.sub("_dataLocAnnot","",ppifiles[i])
ppifiles = list(set(ppifiles))

if not seedRandom == "NA":
	random.seed(seedRandom)

ppifiles = random.sample(ppifiles,permutlength)
print("Building "+str(permutlength)+" permuted networks...")

##make a copy of seeds so it doesn't get altered
seedsCopy = copy.deepcopy(seeds)

for i in range(0,permutlength):
	if not ppifiles[i][0:2]=="RG":
		continue
	permut = buildNet("/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/"+ppifiles[i],"/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/"+ppifiles[i]+"_dataLocAnnot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])
	geneCountDirectPermut = geneCountDirectPermut+(permut["geneCountDirect"],)
	geneCountIndirectPermut = geneCountIndirectPermut+(permut["geneCountIndirect"],)
	locusCountDirectPermut = locusCountDirectPermut+(permut["locusCountDirect"],)
	locusCountIndirectPermut = locusCountIndirectPermut+(permut["locusCountIndirect"],)
	gene2locusEdgeCountDirectPermut = gene2locusEdgeCountDirectPermut+(permut["gene2locusEdgeCountDirect"],)
	gene2locusEdgeCountIndirectPermut = gene2locusEdgeCountIndirectPermut+(permut["gene2locusEdgeCountIndirect"],)
	locus2locusEdgeCountDirectPermut = locus2locusEdgeCountDirectPermut+(permut["locus2locusEdgeCountDirect"],)
	locus2locusEdgeCountIndirectPermut = locus2locusEdgeCountIndirectPermut+(permut["locus2locusEdgeCountIndirect"],)
	gene2locusDegreesMeanDirectPermut = gene2locusDegreesMeanDirectPermut+(permut["gene2locusDegreesMeanDirect"],)
	gene2locusDegreesMeanIndirectPermut = gene2locusDegreesMeanIndirectPermut+(permut["gene2locusDegreesMeanIndirect"],)
	gene2locusWeightedDegreesMeanIndirectPermut = gene2locusWeightedDegreesMeanIndirectPermut+(permut["gene2locusWeightedDegreesMeanIndirect"],)
	locus2locusDegreesMeanDirectPermut = locus2locusDegreesMeanDirectPermut+(permut["locus2locusDegreesMeanDirect"],)
	locus2locusDegreesMeanIndirectPermut = locus2locusDegreesMeanIndirectPermut+(permut["locus2locusDegreesMeanIndirect"],)
	locus2locusWeightedDegreesMeanIndirectPermut = locus2locusWeightedDegreesMeanIndirectPermut+(permut["locus2locusWeightedDegreesMeanIndirect"],)
	CIdegreesMeanPermut = CIdegreesMeanPermut+(permut["CIdegreesMean"],)

	for seed in seeds:
		if permut["gene2locusDegreeDictDirect"].has_key(seed):
			if permut["gene2locusDegreeDictDirect"][seed] >= diseaseNet["gene2locusDegreeDictDirect"][seed]:
				gene2locusDegreeDictDirectPermut[seed]+=1
		if permut["gene2locusDegreeDictIndirect"].has_key(seed):
			if permut["gene2locusDegreeDictIndirect"][seed] >= diseaseNet["gene2locusDegreeDictIndirect"][seed]:
				gene2locusDegreeDictIndirectPermut[seed]+=1
		if permut["gene2locusWeightedDegreeDictIndirect"].has_key(seed):
			if permut["gene2locusWeightedDegreeDictIndirect"][seed] >= diseaseNet["gene2locusWeightedDegreeDictIndirect"][seed]:
				gene2locusWeightedDegreeDictIndirectPermut[seed]+=1
		if permut["locus2locusDegreeDictDirect"].has_key(seed):
			if permut["locus2locusDegreeDictDirect"][seed] >= diseaseNet["locus2locusDegreeDictDirect"][seed]:
				locus2locusDegreeDictDirectPermut[seed]+=1
		if permut["locus2locusDegreeDictIndirect"].has_key(seed):
			if permut["locus2locusDegreeDictIndirect"][seed] >= diseaseNet["locus2locusDegreeDictIndirect"][seed]:
				locus2locusDegreeDictIndirectPermut[seed]+=1
		if permut["locus2locusWeightedDegreeDictIndirect"].has_key(seed):
			if permut["locus2locusWeightedDegreeDictIndirect"][seed] >= diseaseNet["locus2locusWeightedDegreeDictIndirect"][seed]:
				locus2locusWeightedDegreeDictIndirectPermut[seed]+=1

				
	for ci in diseaseNet["CIdegreesDict"].keys():
		if not CIdegreesDictPermut.has_key(ci):
			CIdegreesDictPermut[ci]=0
		if permut["CIdegreesDict"].has_key(ci):
			if permut["CIdegreesDict"][ci] >= diseaseNet["CIdegreesDict"][ci]:
				CIdegreesDictPermut[ci]+=1


	x=range(0,permutlength,10)
	x.append(permutlength-1)
	if i in x:
		print(str(int(float(i)*100/permutlength))+"% complete...")
		tracking=open("tracking"+sys.argv[1],'w')
		tracking.write(str(int(float(i)*99/permutlength))+"% complete...")
		tracking.close()
	    ##write out permuted values
		outTmp = open(str(sys.argv[1])+"_geneCountDirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in geneCountDirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_geneCountIndirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in geneCountIndirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_locusCountDirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in locusCountDirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_gene2locusEdgeCountDirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in gene2locusEdgeCountDirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_gene2locusEdgeCountIndirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in gene2locusEdgeCountIndirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_locus2locusEdgeCountDirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in locus2locusEdgeCountDirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_locus2locusEdgeCountIndirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in locus2locusEdgeCountIndirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_gene2locusDegreesMeanDirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in gene2locusDegreesMeanDirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_gene2locusDegreesMeanIndirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in gene2locusDegreesMeanIndirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_gene2locusWeightedDegreesMeanIndirectPermuted",'w'); outTmp.write("\n".join([str(x) for x in gene2locusWeightedDegreesMeanIndirectPermut])); outTmp.close()
		outTmp = open(str(sys.argv[1])+"_CIdegreesMeanPermuted",'w'); outTmp.write("\n".join([str(x) for x in CIdegreesMeanPermut])); outTmp.close()		


##print result files

print("Summarizing results...")
def P(x,y):
	count = 0
	for num in x:
		if num >= y:
			count+=1
	if count == 0:
		return(float(1)/float(len(x)+1))
	else:
		return((float(count)+1)/(float(len(x))+1))


outfile = open(str(sys.argv[1])+"_GlobalStats",'w')
outfile.write("PARAMETER\tOBSERVED\tEXPECTED\tP_VALUE\n")
outfile.write("geneCountDirect\t"+str(diseaseNet["geneCountDirect"])+"\t"+str(numpy.mean(geneCountDirectPermut))+"\t"+str(P(geneCountDirectPermut,diseaseNet["geneCountDirect"]))+"\n")
outfile.write("geneCountIndirect\t"+str(diseaseNet["geneCountIndirect"])+"\t"+str(numpy.mean(geneCountIndirectPermut))+"\t"+str(P(geneCountIndirectPermut,diseaseNet["geneCountIndirect"]))+"\n")
outfile.write("locusCountDirect\t"+str(diseaseNet["locusCountDirect"])+"\t"+str(numpy.mean(locusCountDirectPermut))+"\t"+str(P(locusCountDirectPermut,diseaseNet["locusCountDirect"]))+"\n")
outfile.write("locusCountIndirect\t"+str(diseaseNet["locusCountIndirect"])+"\t"+str(numpy.mean(locusCountIndirectPermut))+"\t"+str(P(locusCountIndirectPermut,diseaseNet["locusCountIndirect"]))+"\n")
outfile.write("gene2locusEdgeCountDirect\t"+str(diseaseNet["gene2locusEdgeCountDirect"])+"\t"+str(numpy.mean(gene2locusEdgeCountDirectPermut))+"\t"+str(P(gene2locusEdgeCountDirectPermut,diseaseNet["gene2locusEdgeCountDirect"]))+"\n")
outfile.write("gene2locusEdgeCountIndirect\t"+str(diseaseNet["gene2locusEdgeCountIndirect"])+"\t"+str(numpy.mean(gene2locusEdgeCountIndirectPermut))+"\t"+str(P(gene2locusEdgeCountDirectPermut,diseaseNet["gene2locusEdgeCountIndirect"]))+"\n")
outfile.write("locus2locusEdgeCountIndirect\t"+str(diseaseNet["locus2locusEdgeCountIndirect"])+"\t"+str(numpy.mean(locus2locusEdgeCountIndirectPermut))+"\t"+str(P(locus2locusEdgeCountIndirectPermut,diseaseNet["locus2locusEdgeCountIndirect"]))+"\n")
outfile.write("gene2locusDegreesMeanDirect\t"+str(diseaseNet["gene2locusDegreesMeanDirect"])+"\t"+str(numpy.mean(gene2locusDegreesMeanDirectPermut))+"\t"+str(P(gene2locusDegreesMeanDirectPermut,diseaseNet["gene2locusDegreesMeanDirect"]))+"\n")
outfile.write("gene2locusDegreesMeanIndirect\t"+str(diseaseNet["gene2locusDegreesMeanIndirect"])+"\t"+str(numpy.mean(gene2locusDegreesMeanIndirectPermut))+"\t"+str(P(gene2locusDegreesMeanIndirectPermut,diseaseNet["gene2locusDegreesMeanIndirect"]))+"\n")
outfile.write("gene2locusWeightedDegreesMeanIndirect\t"+str(diseaseNet["gene2locusWeightedDegreesMeanIndirect"])+"\t"+str(numpy.mean(gene2locusWeightedDegreesMeanIndirectPermut))+"\t"+str(P(gene2locusWeightedDegreesMeanIndirectPermut,diseaseNet["gene2locusWeightedDegreesMeanIndirect"]))+"\n")
outfile.write("locus2locusDegreesMeanDirect\t"+str(diseaseNet["locus2locusDegreesMeanDirect"])+"\t"+str(numpy.mean(locus2locusDegreesMeanDirectPermut))+"\t"+str(P(locus2locusDegreesMeanDirectPermut,diseaseNet["locus2locusDegreesMeanDirect"]))+"\n")
outfile.write("locus2locusDegreesMeanIndirect\t"+str(diseaseNet["locus2locusDegreesMeanIndirect"])+"\t"+str(numpy.mean(locus2locusDegreesMeanIndirectPermut))+"\t"+str(P(locus2locusDegreesMeanIndirectPermut,diseaseNet["locus2locusDegreesMeanIndirect"]))+"\n")
outfile.write("locus2locusWeightedDegreesMeanIndirect\t"+str(diseaseNet["locus2locusWeightedDegreesMeanIndirect"])+"\t"+str(numpy.mean(locus2locusWeightedDegreesMeanIndirectPermut))+"\t"+str(P(locus2locusWeightedDegreesMeanIndirectPermut,diseaseNet["locus2locusWeightedDegreesMeanIndirect"]))+"\n")
outfile.write("CIdegreesMean\t"+str(diseaseNet["CIdegreesMean"])+"\t"+str(numpy.mean(CIdegreesMeanPermut))+"\t"+str(P(CIdegreesMeanPermut,diseaseNet["CIdegreesMean"]))+"\n")

##Write out seed scores

outfile = open(str(sys.argv[1])+"_seedScores",'w')
outfile.write("SEED\tREGION\tP_gene2locusDirect\tP_gene2locusIndirect\tP_gene2locusIndirectWeighted\n")
countfile = open(str(sys.argv[1])+"_seedPermuteCounts",'w')
##Pdirect is based on seedDirectLocusDegreeDict
##Pindirect is based on seedIndirectWeightedScoreDict
##Assume 2 tests for all seeds
for seed in seedsCopy:
	P_gene2locusDirect = float(gene2locusDegreeDictDirectPermut[seed]+1)/float(permutlength+1)
	P_gene2locusIndirect = float(gene2locusDegreeDictIndirectPermut[seed]+1)/float(permutlength+1)
	P_gene2locusIndirectWeighted = float(gene2locusWeightedDegreeDictIndirectPermut[seed]+1)/float(permutlength+1)
	outfile.write(geneNameConvert[seed]+"\t"+genesnpCopy[seed]+"\t"+str(P_gene2locusDirect)+"\t"+str(P_gene2locusIndirect)+"\t"+str(P_gene2locusIndirectWeighted)+"\n")

outfile.close()

##write out CI scores

outfile = open(sys.argv[1]+"_CIscores",'w')
outfile.write("CI\tNUM_LOCI_CONNECTIONS\tNUM_EXPECTED\tP_VALUE\tP_VALUE_CORRECTED\n")
for ci in diseaseNet["CIdegreesDict"].keys():
	if not ci in CIdegreesDictPermut.keys():
		p = float(1)/(float(permutlength+1)+1)
	else:
		p = (float(CIdegreesDictPermut[ci])+1)/(float(permutlength)+1)
	pfinal = 1-float(1-p)**len(diseaseNet["CIdegreesDict"].keys())
	outfile.write(geneNameConvert[ci]+"\t"+str(len(CIseed[ci]))+"\t"+str(p)+"\t"+str(pfinal)+"\n")

outfile.close()

##write out Genes to Prioritize:

os.system("awk \'$4<0.05{print $1}\' "+sys.argv[1]+"_seedScores > "+sys.argv[1]+"_GenesToPrioritize")

##plot if p-value colored
if plot and plotP:
	nodemap = {}
	nodes = []
	links = []
	links2 = []
	ct = 0
	plotfile = open('data.json','w')
	plotfile2 = open('dataDirect.json','w')
	directSeedList = []
	if plot:
		print("Plotting network...")
		directOut = open("./tmpDirect",'w')
		indirectOut = open("./tmpIndirect",'w')
		for seed in seedSeed:
			source = geneNameConvert[seed]
			if len(seedSeed[seed])>0:
				directSeedList.append(source)
			for otherSeed in seedSeed[seed]:
				pair = [seed,otherSeed]
				pair.sort()
				directOut.write(geneNameConvert[pair[0]]+"\t"+geneNameConvert[pair[1]]+"\n")
				target = geneNameConvert[otherSeed]
				if not nodemap.has_key(source):
					nodemap[source] = ct
					ct += 1
				if not nodemap.has_key(target):
					nodemap[target] = ct
					ct += 1
				d = {
					"source": nodemap[source],
					"target": nodemap[target],
					"value": 2
					}
				links.append(d)
				links2.append(d)
		directOut.close()
		os.system("sort -k 1,1 tmpDirect | uniq > tmp")
		os.system("mv tmp tmpDirect")
		for seed in seedCI:
			source = geneNameConvert[seed]
			for ci in seedCI[seed]:
				target = geneNameConvert[ci]
				pair = [seed,ci]
				pair.sort()
				indirectOut.write(geneNameConvert[pair[0]]+"\t"+geneNameConvert[pair[1]]+"\n")
				if not nodemap.has_key(source):
					nodemap[source] = ct
					ct += 1
				if not nodemap.has_key(target):
					nodemap[target] = ct
					ct += 1					
				d = {
					"source": nodemap[source],
					"target": nodemap[target],
					"value": 1
					}
				links.append(d)
		indirectOut.close()
		os.system("sort -k 1,1 tmpIndirect | uniq > tmp")
		os.system("mv tmp tmpIndirect")
		outfile = open("tmpSnpGene",'w')
		colormap = {}
		col = 1
		snps = list(set(genesnpCopy.values()))
		for seed in genesnpCopy.keys():
			outfile.write(genesnpCopy[seed]+"\t"+geneNameConvert[seed]+"\n")
			colormap[geneNameConvert[seed]] = int(snps.index(genesnpCopy[seed]))
		outfile.close()
		if not collapseCI:
			collapse = "FALSE"
		else:
			collapse = "TRUE"
		if not manualColoring:
			if not plotP:
				os.system("sed s/manualColoring/FALSE/g /web/wwwprod/htdocs/mpg/dapple/NewCode/plot.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")
			else:
				os.system("sed s/manualColoring/FALSE/g /web/wwwprod/htdocs/mpg/dapple/NewCode/plot.R | sed s/NPERMUTEHERE/"+str(permutlength)+"/g | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")			
		else:
			os.system("sed \'s/manualColoring/TRUE/g\' /web/wwwprod/htdocs/mpg/dapple/NewCode/plot.R | sed \'s/coloringScheme/"+manualColoring+"/g\' | sed \'s/KEYWORDHERE/"+sys.argv[1]+"/g\' > tmpPlot.R")
		os.system("R --vanilla < tmpPlot.R")

##get coordinates
	coordFile = open("PlotCoordinates",'r')
	coords = {}
	for line in coordFile:
		line = line.strip("\n").split()
		coords[line[0]] = [float(line[1]),float(line[2])]

##create data for SVG plot
	if len(seeds)<500:
		nodes2=[]
		for k in nodemap.keys():
			color = colormap.get(k, 0)
			d = {
				"index": nodemap[k],
				"name": k,
				"group": color,
				"x": coords[k][0],
				"y": coords[k][1]
				}
			nodes.append(d)
			if k in directSeedList:
				nodes2.append(d)

		nodes = sorted(nodes, key=lambda x: x['index'])
		nodes2 = sorted(nodes2, key=lambda x: x['index'])
		out2 = {
			"nodes": nodes2,
			"links": links2
			}
		json.dump(out2, plotfile2)
		out = {
			"nodes": nodes,
			"links": links
			}
		json.dump(out, plotfile)
		plotfile.close()
		plotfile2.close()

##Iterate, if iterate = 'true'

genestospecify = open(sys.argv[1]+"_GenesToPrioritize",'r').readlines()
for line in genestospecify:
	tokeep.append(geneNameConvert[line.strip("\n")])

tokeep = tokeep+tokeepOriginal
tokeep = list(set(tokeep))
tokeepOriginal = tokeep

seeds=originalseeds[:]
snps=originalsnps[:]


toremove = []
for i in range(len(seeds)):
	if seeds[i] in tokeep:
		toremove.append(snps[i])
indicestoremove = []
for i in range(len(snps)):
	if snps[i] in toremove:
		if seeds[i] in tokeep:
			continue
		else:
			indicestoremove.append(i)

for index in indicestoremove:
	seeds[index] = 'X'
	snps[index] = 'X'
while 'X' in seeds:
	seeds.remove('X')
while 'X' in snps:
	snps.remove('X')

genesnp = {}
### snpGeneCount = {}
for i in range(len(seeds)):
	genesnp[seeds[i]] = snps[i]
###	if not snpGeneCount.has_key(snps[i]):
###		snpGeneCount[snps[i]]=0
###	snpGeneCount[snps[i]]+=1

genesnpCopy = copy.deepcopy(genesnp)
seedsCopy = copy.deepcopy(seeds)
print(seeds)

if not iterate:
	print ("You chose not to iterate. Ending here.")
	sys.exit()

if iterate and len(genestospecify) == 0:
	print ("You chose to iterate but no genes were prioritized for further iteration.")
	sys.exit()

##Iterate until it converges

genesToSpecifyCount = 0
iterateCount = 1
#if iterate and len(genestospecify)>0:
while len(genestospecify) > genesToSpecifyCount:
	genesToSpecifyCount = len(genestospecify)
	iterateCount+=1
	directEdgesCountPermut = ()
	seedDirectDegreesMeanPermut = ()
	seedIndirectDegreesMeanPermut = ()
	CIdegreesMeanPermut = ()
	seedDirectDegreesPermut = {}
	seedIndirectDegreesPermut = {}
	CIdegreesPermut = {}

	ppifiles = os.listdir("/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/")
	for i in range(len(ppifiles)):
		ppifiles[i] = re.sub("_dataLocAnnot","",ppifiles[i])
   	ppifiles = list(set(ppifiles))

	ppifiles = random.sample(ppifiles,permutlength)
	print("Iteration #"+str(iterateCount)+": Building "+str(permutlength)+" permuted networks...")

	print("Building iterated disease network...")

	diseaseNet = buildNet(ppifile,"/home/radon00/rossin/DALY/PPI/NewCode/ppiannot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])
	
	seedSeed = diseaseNet["seedSeed"]
	if any(x!=[] for x in seedSeed.values()):
		count=0
		tally = []
		for seed in seedSeed.keys():
			if len(seedSeed[seed]) > 0:
				tally.append(seed)
				for otherSeed in seedSeed[seed]:
					if not otherSeed in tally:
						count+=1
						print(geneNameConvert[seed]+" -- "+geneNameConvert[otherSeed])
		print("There were "+str(int(float(count)))+" direct interactions in total.")
	else:
		print("There are no direct ineteractions among disease proteins.")

	print("Mean associated protein direct connectivity: "+str(diseaseNet["seedDirectDegreesMean"]))
	print("Mean associated protein indirect connectivity: "+str(diseaseNet["seedIndirectDegreesMean"]))
	print("Mean CI connectivity: "+str(diseaseNet["CIdegreesMean"]))

	for i in range(0,permutlength):
		if not ppifiles[i][0:2]=="RG":
			continue
		permut = buildNet("/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/"+ppifiles[i],"/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/LizzyRandom/ForWebsite/"+ppifiles[i]+"_dataLocAnnot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])
		directEdgesCountPermut = directEdgesCountPermut+(permut["directEdgesCount"],)
		seedDirectDegreesMeanPermut = seedDirectDegreesMeanPermut+(permut["seedDirectDegreesMean"],)
		seedIndirectDegreesMeanPermut = seedIndirectDegreesMeanPermut+(permut["seedIndirectDegreesMean"],)
		CIdegreesMeanPermut = CIdegreesMeanPermut+(permut["CIdegreesMean"],)
		for seed in diseaseNet["seedDirectDegrees"].keys():
			if not seedDirectDegreesPermut.has_key(seed):
				seedDirectDegreesPermut[seed]=0
			if permut["seedDirectDegrees"].has_key(seed):
				if permut["seedDirectDegrees"][seed] >= diseaseNet["seedDirectDegrees"][seed]:
					seedDirectDegreesPermut[seed]+=1
		for seed in diseaseNet["seedIndirectDegrees"].keys():
			if not seedIndirectDegreesPermut.has_key(seed):
				seedIndirectDegreesPermut[seed]=0
			if permut["seedIndirectDegrees"].has_key(seed):
				if permut["seedIndirectDegrees"][seed] >= diseaseNet["seedIndirectDegrees"][seed]:
					seedIndirectDegreesPermut[seed]+=1
		for ci in diseaseNet["CIdegrees"].keys():
			if not CIdegreesPermut.has_key(ci):
				CIdegreesPermut[ci]=0
			if permut["CIdegrees"].has_key(ci):
				if permut["CIdegrees"][ci] >= diseaseNet["CIdegrees"][ci]:
					CIdegreesPermut[ci]+=1
		x=range(0,permutlength,10)
		if i in x:
			print(str(int(float(i)*100/permutlength))+"% complete...")
			tracking=open("tracking",'w')
			tracking.write("Iteration "+str(iterateCount)+" "+str(int(float(i)*100/permutlength))+"% complete...")
			tracking.close()

##print result files

	print("Summarizing iteration #"+str(iterateCount)+" results...")

	outfile = open(str(sys.argv[1])+"_Iterate_NetStats",'w')
	outfile.write("PARAMETER\tP_VALUE\n")
	outfile.write("Direct Edges Count\t"+str(P(directEdgesCountPermut,diseaseNet["directEdgesCount"]))+"\n")
	outfile.write("Seed Direct Degrees Mean\t"+str(P(seedDirectDegreesMeanPermut,diseaseNet["seedDirectDegreesMean"]))+"\n")
	outfile.write("Seed Indirect Degrees Mean\t"+str(P(seedIndirectDegreesMeanPermut,diseaseNet["seedIndirectDegreesMean"]))+"\n")
	outfile.write("CI Degrees Mean\t"+str(P(CIdegreesMeanPermut,diseaseNet["CIdegreesMean"]))+"\n")
	outfile.close()
	
#snpGeneCount is a count of genes per snp

###	snpGeneCount = {}
###	for seed in seeds:
###		if not (diseaseNet["seedDirectDegrees"][seed]==0 and diseaseNet["seedIndirectDegrees"][seed]==0):
###			if not snpGeneCount.has_key(genesnpCopy[seed]):
###				snpGeneCount[genesnpCopy[seed]]=0
###			snpGeneCount[genesnpCopy[seed]]+=1

##Write out seed scores

	outfile = open(str(sys.argv[1])+"_Iterate_seedScores",'w')
	outfile.write("GENE\tREGION\tP_uncorrected\tP_corrected\n")
	for seed in seedsCopy:
		if diseaseNet["seedDirectDegrees"][seed]==0 and diseaseNet["seedIndirectDegrees"][seed]==0:
			outfile.write(geneNameConvert[seed]+"\tNA\tNA\tNA\n")
			continue
		if seed in diseaseNet["seedDirectDegrees"].keys():
			if not seed in seedDirectDegreesPermut.keys():
				pdirect = float(1)/float(permutlength+1)
			else:
				pdirect = (float(seedDirectDegreesPermut[seed])+1)/(float(permutlength)+1)
				if seed in diseaseNet["seedIndirectDegrees"].keys():
					correctby = 2
					if not seed in seedIndirectDegreesPermut.keys():
						pindirect = float(1)/float(permutlength+1)
					else:
						pindirect = (float(seedIndirectDegreesPermut[seed])+1)/(float(permutlength)+1)
					pfinal = min([pdirect,pindirect])
				else:
					correctby = 1
					pfinal = pdirect
		elif seed in diseaseNet["seedIndirectDegrees"].keys():
			correctby = 1
			if not seed in seedIndirectDegreesPermut.keys():
				pindirect = float(1)/float(permutlength+1)
			else:
				pindirect = (float(seedIndirectDegreesPermut[seed])+1)/(float(permutlength)+1)
			pfinal = pindirect
		if pfinal == 0:
			pfinal = float(1)/float(int(permutlength)+1)
		pfinal = 1-float(1-pfinal)**correctby
		correctby=snpGeneCount[genesnpCopy[seed]]
		pfinalCorrected = 1-float(1-pfinal)**correctby
		outfile.write(geneNameConvert[seed]+"\t"+genesnpCopy[seed]+"\t"+str(pfinal)+"\t"+str(pfinalCorrected)+"\n")


	outfile.close()


##write out CI scores

	outfile = open(sys.argv[1]+"_Iterate_CIscores",'w')
	outfile.write("PROTEIN\tNUM_BINDERS\tP_VALUE\tP_VALUE_CORRECTED\n")
	for ci in diseaseNet["CIdegrees"].keys():
		if not ci in CIdegreesPermut.keys():
			p = float(1)/float(permutlength+1)
		else:
			p = float(CIdegreesPermut[ci])/float(permutlength)
			if p==0:
				p = float(1)/float(permutlength+1)
		pfinal = 1-(1-p)**len(diseaseNet["CIdegrees"].keys())
		outfile.write(geneNameConvert[ci]+"\t"+str(len(CIseed[ci]))+"\t"+str(p)+"\t"+str(pfinal)+"\n")

	outfile.close()

##write out Genes to Prioritize:

	os.system("awk \'$4<0.05{print $1}\' "+sys.argv[1]+"_Iterate_seedScores > "+sys.argv[1]+"_Iterate_GenesToPrioritize")

##organize any new genes that are prioritized
	genestospecify = open(sys.argv[1]+"_Iterate_GenesToPrioritize",'r').readlines()
	for line in genestospecify:
		tokeep.append(geneNameConvert[line.strip("\n")])

	tokeep = tokeep+tokeepOriginal
	tokeep = list(set(tokeep))
	tokeepOriginal = tokeep
	
	seeds=originalseeds[:]
	snps=originalsnps[:]

	toremove = []
	for i in range(len(seeds)):
		if seeds[i] in tokeep:
			toremove.append(snps[i])
	indicestoremove = []
	for i in range(len(snps)):
		if snps[i] in toremove:
			if seeds[i] in tokeep:
				continue
			else:
				indicestoremove.append(i)

	for index in indicestoremove:
		seeds[index] = 'X'
		snps[index] = 'X'
	while 'X' in seeds:
		seeds.remove('X')
	while 'X' in snps:
		snps.remove('X')

	genesnp = {}
###	snpGeneCount = {}
	for i in range(len(seeds)):
		genesnp[seeds[i]] = snps[i]
###		if not snpGeneCount.has_key(snps[i]):
###			snpGeneCount[snps[i]]=0
###		snpGeneCount[snps[i]]+=1

	genesnpCopy = copy.deepcopy(genesnp)
	seedsCopy = copy.deepcopy(seeds)
	
print("Iteration converged. Exiting.")
