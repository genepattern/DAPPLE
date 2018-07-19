
import os
import sys
import re
import random
import copy
import json
import numpy
from Functions_STRING import *

##input can be SNPs, genes or regions

inputlist = open(sys.argv[1],'r').readlines()
genestospecify = "X"
CIcutoff = 2
permutlength=1000
plot = 'false'
plotP = 'false'
manualColoring = 'false'
iterate = 'false'
seedprint = 'false'
compare = 'false'
regulatory = [50000,50000]
dir = "./"
nearestgene = "false"
seedRandom = "NA"
collapseCI = "false"
plotP = "false"
if len(sys.argv)>2:
	for argv in sys.argv[2:]:
		if argv.split("=")[0] == "specify":
			genestospecify = open(argv.split("=")[1],'r')
		elif argv.split("=")[0] == "permute":
			permutlength = int(argv.split("=")[1])
		elif argv.split("=")[0] == "CIcutoff":
			CIcutoff = int(argv.split("=")[1])
		elif argv.split("=")[0] == "plot":
			plot = argv.split("=")[1]
		elif argv.split("=")[0] == "iterate":
			iterate = argv.split("=")[1]
		elif argv.split("=")[0] == "seedprint":
			seedprint = argv.split("=")[1]
		elif argv.split("=")[0] == "manual":
			manualColoring = argv.split("=")[1]
		elif argv.split("=")[0] == "keyword":
			sys.argv[1] = argv.split("=")[1]
		elif argv.split("=")[0] == "compare":
			compare = argv.split("=")[1]
		elif argv.split("=")[0] == "regulatory":
			regulatory = [1000*int(argv.split("=")[1].split(",")[0]),1000*int(argv.split("=")[1].split(",")[1])]
		elif argv.split("=")[0] == "dir":
			dir = argv.split("=")[1]
		elif argv.split("=")[0] == "nearestgene":
			nearestgene = argv.split("=")[1]
		elif argv.split("=")[0] == "seed":
			seedRandom = int(argv.split("=")[1])
		elif argv.split("=")[0] == "plotP":
			plotP = argv.split("=")[1]
		elif argv.split("=")[0] == "collapseCI":
			collapseCI = argv.split("=")[1]
		elif argv.split("=")[0] == "database":
			database = argv.split("=")[1]

##Change to the specified working directory
os.chdir(dir)

##load PPI data and create dictionaries
if database == "full":
	ppifile = "/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/STRING_long"
if database == "400":
	ppifile = "/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/STRING_heiko_400_underscore"
if database == "700":
	ppifile = "/fg/wgas2/rossin/KasperData/Diana/RandomGraphs/STRING_heiko_700_underscore"

##create a list of proteins in the high-confidence set
hcProteins = []
ppi = open(ppifile,'r')
for line in ppi:
	line = line.strip("\n").split("-")
	hcProteins.append(line[0])
	hcProteins.append(line[1])
hcProteins = list(set(hcProteins))

##set gene input as false
geneInput = False

##if input is a list of SNPs and nearestgene = true
if nearestgene == "true":
	print("Input is a list of snps and you have chosen to find the nearest gene")
	snps = []
	for line in inputlist:
		snps.append(line.strip("\n").strip("\r"))
	x=NearestGene(snps)
	seeds = x[0]
	snps = x[1]
##if input is a list of SNPs
elif inputlist[0][0:2]=="rs" and len(inputlist[0].split())==1:
	print("Input is a list of snps. Converting snps to genes...")
	snps = []
	for line in inputlist:
		snps.append(line.strip("\n").strip("\r"))
	if compare == 'false':
		x = SNPtoGene(snps,regulatory[0],regulatory[1])
	elif compare == 'true':
		x = SNPtoGeneCompare(snps,regulatory[0],regulatory[1])
	seeds = x[0]
	snps = x[1]

##if input is a list of genes
elif (len(inputlist[0].split())==1 and inputlist[0][0:2]!="rs") or len(inputlist[0].split())==2:
	print("Input is a list of genes. Extracting genes...")
	inputlist = list(set(inputlist))
	geneInput = True
	count=0
	seeds=[]
	snps=[]
	missingSeeds=[]
	missingSNPs=[]
	for line in inputlist:
		if line.strip("\n").split()[0] in hcProteins:
			seeds.append(line.strip("\n").split()[0])
			if len(inputlist[0].split())==2:
				snps.append(line.strip("\n").split()[1])
			else:
				snps.append(str(count))
				count+=1
		else:
			missingSeeds.append(line.strip("\n").split()[0])
			if len(inputlist[0].split())==2:
				missingSNPs.append(line.strip("\n").split()[1])
			else:
				missingSNPs.append(str(count))
				count+=1

##if input is a list of regions
elif len(inputlist[0].split())==4:
	print("Input is a list of regions. Converting regions to genes...")
	x = RegionToGene(inputlist,regulatory[0],regulatory[1])
	seeds = x[0]
	snps = x[1]
else:
	print("Input format is invalid. Input must be one of the following:\n-SNP list (one SNP per line, each identified by an rs number)\n-Gene list (each line consisting of one gene or one gene and a region ID, separated by a tab)\n-One region per line (each line consists of a region ID, chromosome, left bound and right bound, all tab separated)")
	sys.exit()

##Exit if too few genes match PPI data
if float(len(seeds))/float(len(inputlist)) < 0.1:
	print("Less tahn 10% of seed genes are in the PPI database. That is too few. Check to make sure IDs are Hugo IDs - for example, \"ATXN1\". If your input was a list of SNPs, make sure they are in HapMap release 22.")
	sys.exit()
	
##Exit if seedprint = 'true'
if seedprint == 'true':
	print("Seeds are:")
	for i in range(len(seeds)):
		print(seeds[i]+"\t"+snps[i])
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
   		tokeep.append(line.strip("\n"))
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

##write out missing genes file 
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
	outfile = open(sys.argv[1]+"_MissingGenes",'w')
	outfile.write("SNP\tALLCount\tPPICount\tPercent\tPPIGenes\tMissingGenes\n")
	allSnps = list(set(snps+missingSNPs))
	for snp in allSnps:
		outfile.write(snp+"\t")
		outfile.write(str(snpGeneCount[snp]+missingGeneCount[snp])+"\t"+str(snpGeneCount[snp])+"\t")
		outfile.write(str(float(snpGeneCount[snp])/float(snpGeneCount[snp]+missingGeneCount[snp]))+"\t")
		for gene in genesnp.keys():
			if genesnp[gene]==snp:
				outfile.write(gene+",")
		outfile.write("\t")
		for gene in missingGeneSnp.keys():
			if missingGeneSnp[gene]==snp:
				outfile.write(gene+",")
		outfile.write("\n")

##print seed-locus file
seedLocusFile = open(sys.argv[1]+"_seedLocusMapping",'w')
for seed in seeds:
	seedLocusFile.write(seed+"\t"+genesnp[seed]+"\n")
seedLocusFile.close()


##run disease network, set up dictionaries
print("Building disease network...")


diseaseNet = buildNet(ppifile,ppifile+"_dataLocAnnot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])

seedSeed = diseaseNet["seedSeed"]
seedCI = diseaseNet["seedCI"]
CIseed = diseaseNet["CIseed"]

outfile = open(sys.argv[1]+"_directConnections",'w')
count=0
if any(x!=[] for x in seedSeed.values()):
	count=0
	tally = []
	for seed in seedSeed.keys():
		if len(seedSeed[seed]) > 0:
			tally.append(seed)
			for otherSeed in seedSeed[seed]:
				if not otherSeed in tally:
					count+=1
					print(seed+" -- "+otherSeed)
					outfile.write(seed+" -- "+otherSeed+"\n")
	print("There were "+str(int(float(count)))+" direct interactions in total.")
else:
	print("There are no direct ineteractions among disease proteins.")
outfile.close()
outfile = open(sys.argv[1]+"_summary",'w')
if diseaseNet["seedDirectDegreesMean"]==0:
	outfile.write("There were "+str(0)+" disease proteins participating in the direct network.\n")	
	print("Number of disease proteins in direct network: "+str(0))
else:
	outfile.write("There were "+str(int(float(count))*2/diseaseNet["seedDirectDegreesMean"])+" disease proteins participating in the direct network.\n")
	print("Number of disease proteins in direct network: "+str(int(float(count))*2/diseaseNet["seedDirectDegreesMean"]))
outfile.write("There were "+str(int(float(count)))+" direct interactions in total."+"\n")
print("Mean associated protein direct connectivity: "+str(diseaseNet["seedDirectDegreesMean"]))
outfile.write("Mean associated protein direct connectivity: "+str(diseaseNet["seedDirectDegreesMean"])+"\n")
print("Mean associated protein indirect connectivity: "+str(diseaseNet["seedIndirectDegreesMean"]))
outfile.write("Mean associated protein indirect connectivity: "+str(diseaseNet["seedIndirectDegreesMean"])+"\n")
print("Mean CI connectivity: "+str(diseaseNet["CIdegreesMean"])) 
outfile.write("Mean CI connectivity: "+str(diseaseNet["CIdegreesMean"])+"\n")
outfile.close()

if count == 0 and diseaseNet["seedIndirectDegreesMean"] == 0:
	plot='false'

	
outfile = open(sys.argv[1]+"_indirectConnections",'w')
for seed in seedCI.keys():
	if (len(seedCI[seed])>0):
		for CI in seedCI[seed]:
			outfile.write(seed+" -- "+CI+"\n")
outfile.close()


if (plot == 'true' and plotP == 'false') or (plot == 'true' and permutlength==0):
	nodemap = {}
	nodes = []
	links = []
	links2 = []
	ct = 0
	plotfile = open('data.json','w')
	plotfile2 = open('dataDirect.json','w')
	directSeedList = []
	if plot=='true':
		print("Plotting network...")
		directOut = open("./tmpDirect",'w')
		indirectOut = open("./tmpIndirect",'w')
		for seed in seedSeed:
			source = seed
			if len(seedSeed[seed])>0:
				directSeedList.append(source)
			for otherSeed in seedSeed[seed]:
				pair = [seed,otherSeed]
				pair.sort()
				directOut.write(pair[0]+"\t"+pair[1]+"\n")
				target = otherSeed
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
			source = seed
			for ci in seedCI[seed]:
				target = ci
				pair = [seed,ci]
				pair.sort()
				indirectOut.write(pair[0]+"\t"+pair[1]+"\n")
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
			outfile.write(genesnp[seed]+"\t"+seed+"\n")
			colormap[seed] = snps.index(genesnp[seed])
		outfile.close()
		if collapseCI == 'false':
			collapse = "FALSE"
		else:
			collapse = "TRUE"
		if manualColoring == 'false':
			if plotP == 'false':
				os.system("sed s/manualColoring/FALSE/g /home/radon00/rossin/DALY/PPI/NewCode/plot.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")
			else:
				os.system("sed s/manualColoring/FALSE/g /home/radon00/rossin/DALY/PPI/NewCode/plotP.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")			
		else:
			os.system("sed \'s/manualColoring/TRUE/g\' /home/radon00/rossin/DALY/PPI/NewCode/plot.R | sed \'s/coloringScheme/"+manualColoring+"/g\' | sed \'s/KEYWORDHERE/"+sys.argv[1]+"/g\' > tmpPlot.R")
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


##exit if permutlength=0
if permutlength==0:
	sys.exit()

##permute
directEdgesCountPermut = ()
seedDirectDegreesMeanPermut = ()
seedIndirectDegreesMeanPermut = ()
CIdegreesMeanPermut = ()
seedDirectDegreesPermut = {}
seedIndirectDegreesPermut = {}
CIdegreesPermut = {}

##select database
if database == "full":
	directory = "/humgen/atgu1/fs03/wip/rossin/STRINGrandom_heiko/"
elif database == "400":
	directory = "/humgen/atgu1/fs03/wip/rossin/STRINGrandom_heiko_400/"
elif database == "700":
	directory = "/humgen/atgu1/fs03/wip/rossin/STRINGrandom_heiko_700/"

ppifiles = os.listdir(directory)
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
	permut = buildNet(directory+ppifiles[i],directory+ppifiles[i]+"_dataLocAnnot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])
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
		tracking.write("Iteration 1 "+str(int(float(i)*100/permutlength))+"% complete...")
		tracking.close()

##write out permuted values
outTmp = open(str(sys.argv[1])+"_permutedDirectEdgeCount",'w')
for count in directEdgesCountPermut:
	outTmp.write(str(count)+"\n")
outTmp.close()
outTmp = open(str(sys.argv[1])+"_permutedSeedDirectDegreesMean",'w')
for count in seedDirectDegreesMeanPermut:
	outTmp.write(str(count)+"\n")
outTmp.close()
outTmp = open(str(sys.argv[1])+"_permutedIndirectDegreesMean",'w')
for count in seedIndirectDegreesMeanPermut:
	outTmp.write(str(count)+"\n")
outTmp.close()
outTmp = open(str(sys.argv[1])+"_permutedCIdegreesMeanPermut",'w')
for count in CIdegreesMeanPermut:
	outTmp.write(str(count)+"\n")
outTmp.close()
outTmp = open(str(sys.argv[1])+"_permutedNumDirectProteins",'w')
for i in range(0,len(directEdgesCountPermut)):
	if seedDirectDegreesMeanPermut[i]==0:
		outTmp.write("0\n")
	else:
		outTmp.write(str(directEdgesCountPermut[i]*2/seedDirectDegreesMeanPermut[i])+"\n")
outTmp.close()



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

outfile = open(str(sys.argv[1])+"_NetStats",'w')
outfile.write("PARAMETER\tOBSERVED\tEXPECTED\tP_VALUE\n")
outfile.write("Direct Edges Count\t"+str(diseaseNet["directEdgesCount"])+"\t"+str(numpy.mean(directEdgesCountPermut))+"\t"+str(P(directEdgesCountPermut,diseaseNet["directEdgesCount"]))+"\n")
outfile.write("Seed Direct Degrees Mean\t"+str(diseaseNet["seedDirectDegreesMean"])+"\t"+str(numpy.mean(seedDirectDegreesMeanPermut))+"\t"+str(P(seedDirectDegreesMeanPermut,diseaseNet["seedDirectDegreesMean"]))+"\n")
outfile.write("Seed Indirect Degrees Mean\t"+str(diseaseNet["seedIndirectDegreesMean"])+"\t"+str(numpy.mean(seedIndirectDegreesMeanPermut))+"\t"+str(P(seedIndirectDegreesMeanPermut,diseaseNet["seedIndirectDegreesMean"]))+"\n")
outfile.write("CI Degrees Mean\t"+str(diseaseNet["CIdegreesMean"])+"\t"+str(numpy.mean(CIdegreesMeanPermut))+"\t"+str(P(CIdegreesMeanPermut,diseaseNet["CIdegreesMean"]))+"\n")
outfile.close()

#write out direct connections
#directOutCount = open("permutedEdgeCount",'w')
#for num in directEdgesCountPermut:
#	directOutCount.write(str(num)+"\n")
#directOutCount.close()

#snpGeneCount is a count of genes per snp

###snpGeneCount = {}
###for seed in seedsCopy:
###	if not (diseaseNet["seedDirectDegrees"][seed]==0 and diseaseNet["seedIndirectDegrees"][seed]==0):
###		if not snpGeneCount.has_key(genesnpCopy[seed]):
###			snpGeneCount[genesnpCopy[seed]]=0
###		snpGeneCount[genesnpCopy[seed]]+=1

##Write out seed scores

outfile = open(str(sys.argv[1])+"_SeedScores",'w')
outfile.write("GENE\tREGION\tP_uncorrected\tP_corrected\n")

for seed in seedsCopy:
	if diseaseNet["seedDirectDegrees"][seed]==0 and diseaseNet["seedIndirectDegrees"][seed]==0:
		outfile.write(seed+"\tNA\tNA\tNA\n")
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
	outfile.write(seed+"\t"+genesnpCopy[seed]+"\t"+str(pfinal)+"\t"+str(pfinalCorrected)+"\n")

outfile.close()


##write out CI scores

outfile = open(sys.argv[1]+"_CIscores",'w')
outfile.write("PROTEIN\tNUM_BINDERS\tP_VALUE\tP_VALUE_CORRECTED\n")
for ci in diseaseNet["CIdegrees"].keys():
	if not ci in CIdegreesPermut.keys():
		p = float(1)/float(permutlength+1)
	else:
		p = float(CIdegreesPermut[ci])/float(permutlength)
		if p == 0:
			p = float(1)/float(permutlength+1)
	pfinal = 1-float(1-p)**len(diseaseNet["CIdegrees"].keys())
	outfile.write(ci+"\t"+str(len(CIseed[ci]))+"\t"+str(p)+"\t"+str(pfinal)+"\n")

outfile.close()

##write out Genes to Prioritize:

os.system("awk \'$4<0.05{print $1}\' "+sys.argv[1]+"_SeedScores > "+sys.argv[1]+"_GenesToPrioritize")

##plot if p-value colored
if plot == 'true' and plotP == 'true':
	nodemap = {}
	nodes = []
	links = []
	links2 = []
	ct = 0
	plotfile = open('data.json','w')
	plotfile2 = open('dataDirect.json','w')
	directSeedList = []
	if plot=='true':
		print("Plotting network...")
		directOut = open("./tmpDirect",'w')
		indirectOut = open("./tmpIndirect",'w')
		for seed in seedSeed:
			source = seed
			if len(seedSeed[seed])>0:
				directSeedList.append(source)
			for otherSeed in seedSeed[seed]:
				pair = [seed,otherSeed]
				pair.sort()
				directOut.write(pair[0]+"\t"+pair[1]+"\n")
				target = otherSeed
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
			source = seed
			for ci in seedCI[seed]:
				target = ci
				pair = [seed,ci]
				pair.sort()
				indirectOut.write(pair[0]+"\t"+pair[1]+"\n")
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
			outfile.write(genesnpCopy[seed]+"\t"+seed+"\n")
			colormap[seed] = snps.index(genesnpCopy[seed])
		outfile.close()
		if collapseCI == 'false':
			collapse = "FALSE"
		else:
			collapse = "TRUE"
		if manualColoring == 'false':
			if plotP == 'false':
				os.system("sed s/manualColoring/FALSE/g /home/radon00/rossin/DALY/PPI/NewCode/plot.R | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")
			else:
				os.system("sed s/manualColoring/FALSE/g /home/radon00/rossin/DALY/PPI/NewCode/plotP.R | sed s/NPERMUTEHERE/"+str(permutlength)+"/g | sed s/KEYWORDHERE/"+sys.argv[1]+"/g | sed s/COLLAPSEHERE/"+collapse+"/g > tmpPlot.R")			
		else:
			os.system("sed \'s/manualColoring/TRUE/g\' /home/radon00/rossin/DALY/PPI/NewCode/plot.R | sed \'s/coloringScheme/"+manualColoring+"/g\' | sed \'s/KEYWORDHERE/"+sys.argv[1]+"/g\' > tmpPlot.R")
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
	tokeep.append(line.strip("\n"))

tokeep = tokeep+tokeepOriginal
tokeep = list(set(tokeep))
tokeepOriginal = tokeep

seeds=originalseeds[:]
snps=originalsnps[:]
print(seeds)

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
###snpGeneCount = {}
for i in range(len(seeds)):
	genesnp[seeds[i]] = snps[i]
###	if not snpGeneCount.has_key(snps[i]):
###		snpGeneCount[snps[i]]=0
###	snpGeneCount[snps[i]]+=1

genesnpCopy = copy.deepcopy(genesnp)
seedsCopy = copy.deepcopy(seeds)
print(seeds)

if iterate == 'false':
	print ("You chose not to iterate. Ending here.")
	sys.exit()

if iterate == 'true' and len(genestospecify) == 0:
	print ("You chose to iterate but no genes were prioritized for further iteration.")
	sys.exit()

##Iterate until it converges

genesToSpecifyCount = 0
iterateCount = 1
#if iterate == 'true' and len(genestospecify)>0:
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

	ppifiles = os.listdir("/humgen/atgu1/fs03/wip/rossin/STRINGrandom_heiko/")
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
						print(seed+" -- "+otherSeed)
		print("There were "+str(int(float(count)))+" direct interactions in total.")
	else:
		print("There are no direct ineteractions among disease proteins.")

	print("Mean associated protein direct connectivity: "+str(diseaseNet["seedDirectDegreesMean"]))
	print("Mean associated protein indirect connectivity: "+str(diseaseNet["seedIndirectDegreesMean"]))
	print("Mean CI connectivity: "+str(diseaseNet["CIdegreesMean"]))

	for i in range(0,permutlength):
		if not ppifiles[i][0:2]=="RG":
			continue
		permut = buildNet(directory+ppifiles[i],directory+ppifiles[i]+"_dataLocAnnot",seeds,genesnp,CIcutoff,keyword=sys.argv[1])
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

	outfile = open(str(sys.argv[1])+"_Iterate_SeedScores",'w')
	outfile.write("GENE\tREGION\tP_uncorrected\tP_corrected\n")
	for seed in seedsCopy:
		if diseaseNet["seedDirectDegrees"][seed]==0 and diseaseNet["seedIndirectDegrees"][seed]==0:
			outfile.write(seed+"\tNA\tNA\tNA\n")
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
		outfile.write(seed+"\t"+genesnpCopy[seed]+"\t"+str(pfinal)+"\t"+str(pfinalCorrected)+"\n")


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
		outfile.write(ci+"\t"+str(len(CIseed[ci]))+"\t"+str(p)+"\t"+str(pfinal)+"\n")

	outfile.close()

##write out Genes to Prioritize:

	os.system("awk \'$4<0.05{print $1}\' "+sys.argv[1]+"_Iterate_SeedScores > "+sys.argv[1]+"_Iterate_GenesToPrioritize")

##organize any new genes that are prioritized
	genestospecify = open(sys.argv[1]+"_Iterate_GenesToPrioritize",'r').readlines()
	for line in genestospecify:
		tokeep.append(line.strip("\n"))

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
