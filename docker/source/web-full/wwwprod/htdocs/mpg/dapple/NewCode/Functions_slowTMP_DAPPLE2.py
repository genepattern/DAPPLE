import itertools

def buildNet(ppifile,dataAnnotFile,seeds,genesnp,CIcutoff,keyword):
	import copy
	ppifile = open(ppifile,'r').readlines()
	ppiannot = open(dataAnnotFile,'r').readlines()
	dataLoc = {}
	for line in ppiannot:
		line = line.strip("\n").split("\t")
		dataLoc[int(line[0])] = (int(line[1]),int(line[2]))
	seedFind = {}
	for line in ppiannot:
		line = line.strip("\n").split("\t")
		seedFind[int(line[0])]=int(line[1])

	seedSeedL = {}
	seedCI = {}
	seedSeed = {}
	seedCIseedL = {}
	
	for seed in seeds:
		seedCI[seed] = []
		seedSeedL[seed] = []	
		seedSeed[seed] = []
		seedCIseedL[seed] = []

	CIseed = {}
	CIseedL = {}
	

	if "\t" in ppifile[1]:
		d = "\t"
	elif "-" in ppifile[1]:
		d = "-"
	else:
		print(ppifile[1])
		print(" does not have a delimiter")

##begin time consuming part (rate-limiting)

	seedscopy = copy.deepcopy(seeds)
	for seed in seedscopy:
		if not dataLoc.has_key(seed):
			seeds.remove(seed)			
			del genesnp[seed]
			continue
		##Get submatrix with first column as seed
		subData = ppifile[dataLoc[seed][0]:dataLoc[seed][1]]
		for line in subData:
			string_pair = line.strip("\n").split(d)
			pair = (int(string_pair[0]), int(string_pair[1]))
			##Look for second column seeds and add to seedSeed
			if pair[1] in seeds:
				if genesnp[pair[0]] != genesnp[pair[1]]:
					seedSeed[pair[0]].append(pair[1])
					if not genesnp[pair[1]] in seedSeedL[pair[0]]:
						seedSeedL[pair[0]].append(genesnp[pair[1]])
			##Second column always acts as a CI
			##Only evaluate if it hasn't been seen before
			##CI has to have connections to >1 loci
			CI = pair[1]
			if CI in CIseed:
				continue
			subData2 = ppifile[dataLoc[CI][0]:dataLoc[CI][1]]
			tmpCIseedList = []
			tmpCIseedLList = []
			for line2 in subData2:
				line2 = line2.strip("\n").split(d)
				pair2 = (int(line2[0]), int(line2[1]))
				if pair2[1] in seeds and pair2[1] != pair[0]:
					tmpCIseedList.append(pair2[1])
					tmpCIseedLList.append(genesnp[pair2[1]])
					if genesnp[seed] != genesnp[pair2[1]]:
						seedCIseedL[seed].append(genesnp[pair2[1]]+"-"+str(seed))
			if len(list(set(tmpCIseedLList))) > 1:
				CIseed[CI] = tmpCIseedList
				CIseedL[CI] = list(set(tmpCIseedLList))
				seedCI[pair[0]].append(CI)



##end time consuming part

##Count the number of subgraphs in direct network

#	subGraphs = []
#	tally=[]
#	tmpnet=[]
#	for seed in seedSeed.keys():
#		if not seed in tally and len(seedSeed[seed])>0:
#			tmpnet = [seed]+seedSeed[seed]
#			size = len(tmpnet)
#            size2 = 0
#			while size2!=size:
#				size = len(tmpnet)
#				tmpnetcopy = tmpnet[:]
#				for otherSeed in tmpnetcopy:
#					tmpnet = tmpnet+seedSeed[otherSeed]
#					tmpnet = list(set(tmpnet))
#					tally = tally+tmpnet
#				size2 = len(tmpnet)
#			subGraphs.append(tmpnet)

##Create degrees dictionaries, count statistics

	geneCountDirect = 0
	geneCountIndirect = 0
	locusListDirect = []; locusCountDirect = 0
	locusListIndirect = []; locusCountIndirect = 0
	gene2locusEdgeCountDirect = 0
	gene2locusEdgeCountIndirect = 0
	locus2locusEdgeListDirect = []
	locus2locusEdgeCountDirect = 0
	locus2locusEdgeListIndirect = []
	locus2locusEdgeCountIndirect = 0
	gene2locusDegreeDictDirect = {}
	locus2locusDegreeDictDirect = {}
	gene2locusDegreeDictIndirect = {}
	locus2locusDegreeDictIndirect = {}
	gene2locusWeightedDegreeDictIndirect = {}
	locus2locusWeightedDegreeDictIndirect = {}
	
	directTally=[]
	indirectTally = []

	for seed in seeds:
		gene2locusDegreeDictDirect[seed] = len(seedSeedL[seed])
		gene2locusEdgeCountDirect += len(seedSeedL[seed])
		if len(seedCI[seed]) > 0:
			geneCountIndirect+=1
			locusListIndirect.append(genesnp[seed])
			loci = []
			for ci in seedCI[seed]:
				loci = loci+CIseedL[ci]
				x=CIseedL[ci]
				x.sort()
				locus2locusEdgeListIndirect = locus2locusEdgeListIndirect+list(itertools.combinations(x,2))
			gene2locusDegreeDictIndirect[seed] = len(list(set(loci)))
			gene2locusEdgeCountIndirect += len(list(set(loci)))
		else:
			gene2locusDegreeDictIndirect[seed] = 0
		if len(seedSeed[seed]) > 0:
			geneCountDirect+=1
			locusListDirect.append(genesnp[seed])
			for otherLocus in seedSeedL[seed]:
				y = (genesnp[seed],otherLocus)
				y = tuple(sorted(y))
				locus2locusEdgeListDirect.append(y)

	locus2locusEdgeListDirect = list(set(locus2locusEdgeListDirect))
	locus2locusEdgeListIndirect = list(set(locus2locusEdgeListIndirect))
	for edge in locus2locusEdgeListDirect:
		if not locus2locusDegreeDictDirect.has_key(edge[0]):
			locus2locusDegreeDictDirect[edge[0]] = 0
		if not locus2locusDegreeDictDirect.has_key(edge[1]):
			locus2locusDegreeDictDirect[edge[1]] = 0
		locus2locusDegreeDictDirect[edge[0]]+=1
		locus2locusDegreeDictDirect[edge[1]]+=1
	for edge in locus2locusEdgeListIndirect:
		if not locus2locusDegreeDictIndirect.has_key(edge[0]):
			locus2locusDegreeDictIndirect[edge[0]] = 0
		if not locus2locusDegreeDictIndirect.has_key(edge[1]):
			locus2locusDegreeDictIndirect[edge[1]] = 0
		locus2locusDegreeDictIndirect[edge[0]]+=1
		locus2locusDegreeDictIndirect[edge[1]]+=1
	
	allLoci = list(set(genesnp.values()))
	for seed in seeds:
		gene2locusWeightedDegreeDictIndirect[seed] = 0
		if not locus2locusWeightedDegreeDictIndirect.has_key(genesnp[seed]):
			locus2locusWeightedDegreeDictIndirect[genesnp[seed]] = 0
		tmpLocusVec = seedCIseedL[seed]
		for locus in allLoci:
			subset = [i for i in tmpLocusVec if locus in i]
			indirectConnectValues = [0]
			for target in list(set(subset)):
				indirectConnectValues.append(subset.count(target))
			gene2locusWeightedDegreeDictIndirect[seed] += 1-0.9**max(indirectConnectValues)
		if gene2locusWeightedDegreeDictIndirect[seed] > locus2locusWeightedDegreeDictIndirect[genesnp[seed]]:
			locus2locusWeightedDegreeDictIndirect[genesnp[seed]] = gene2locusWeightedDegreeDictIndirect[seed]
		

	CIdegreesDict = {}
	for ci in CIseed.keys():
		CIdegreesDict[ci]=len(list(set(CIseedL[ci])))

##summarize results
	def mean(x):
		if len(x)==0:
			return(0)
		return(float(sum(x))/float(len(x)))

##return results

	results={}
	results["geneCountDirect"] = geneCountDirect
	results["geneCountIndirect"] = geneCountIndirect
	results["locusCountDirect"] = len(list(set(locusListDirect)))
	results["locusCountIndirect"] = len(list(set(locusListIndirect)))
	results["gene2locusEdgeCountDirect"] = gene2locusEdgeCountDirect
	results["gene2locusEdgeCountIndirect"] = gene2locusEdgeCountIndirect
	results["locus2locusEdgeCountDirect"] = len(locus2locusEdgeListDirect)
	results["locus2locusEdgeCountIndirect"] = len(locus2locusEdgeListIndirect)
	results["gene2locusDegreeDictDirect"] = gene2locusDegreeDictDirect
	results["gene2locusDegreesMeanDirect"] = mean(gene2locusDegreeDictDirect.values())
	results["gene2locusDegreeDictIndirect"] = gene2locusDegreeDictIndirect
	results["gene2locusDegreesMeanIndirect"] = mean(gene2locusDegreeDictIndirect.values())
	results["gene2locusWeightedDegreeDictIndirect"] = gene2locusWeightedDegreeDictIndirect
	results["gene2locusWeightedDegreesMeanIndirect"] = mean(gene2locusWeightedDegreeDictIndirect.values())
	results["locus2locusDegreeDictDirect"] = locus2locusDegreeDictDirect
	results["locus2locusDegreesMeanDirect"] = mean(locus2locusDegreeDictDirect.values())
	results["locus2locusDegreeDictIndirect"] = locus2locusDegreeDictIndirect
	results["locus2locusDegreesMeanIndirect"] = mean(locus2locusDegreeDictIndirect.values())
	results["locus2locusWeightedDegreeDictIndirect"] = locus2locusWeightedDegreeDictIndirect
	results["locus2locusWeightedDegreesMeanIndirect"] = mean(locus2locusWeightedDegreeDictIndirect.values())
   	results["CIdegreesDict"] = CIdegreesDict
	results["CIdegreesMean"] = mean(CIdegreesDict.values())
	results["seedSeed"] = seedSeed
	results["seedCI"] = seedCI
	results["CIseed"] = CIseed
#	results["subGraphs"] = subGraphs
	return(results)



#####SNP to Gene Function#######

###Define list of proteins in the InWeb3 HC set
HC=[]
HCfile = open("/fg/wgas2/rossin/KasperData/InWeb3_HC_Red_hugo",'r').readlines()
for line in HCfile:
	HC.append(line.strip("\n").split()[0])
	HC.append(line.strip("\n").split()[1])
HC = list(set(HC))


def SNPtoGeneOLD(snps,regUp,regDown,build):
	import os
	import sys
	import re
	if build == "18":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg18.txt",'r').readlines()

	if build == "19" or build == "1kg":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg19.txt",'r').readlines()
	genes = []
	genesnps = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		if build == "18":
			wsfile = open("/fg/wgas2/rossin/split/SNAPLDFiles/chr"+str(i)+".ws",'r').readlines()
		if build == "19" or build == "1kg":
			wsfile = open("/fg/wgas2/rossin/split/hg19/chr"+str(i)+".ws",'r').readlines()
		for line in wsfile:
			if line.split("\t")[0] in snps:
				count+=1
				snp = line.split("\t")[0]
				chr = str(i)
				left = float(line.split("\t")[1])
				right = float(line.strip("\n").split("\t")[2])
				for geneline in genelist:
					geneline = geneline.strip("\n").split("\t")
					if geneline[2]==chr:
						if geneline[3]=="+":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								genes.append(int(re.sub(".p","",geneline[12])))
								genesnps.append(line.split("\t")[0])
						if geneline[3]=="-":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								genes.append(int(re.sub(".p","",geneline[12])))
								genesnps.append(line.split("\t")[0])
	results={}
	results[0]=genes
	results[1]=genesnps
	return(results)


def SNPtoGene(snps,regUp,regDown,build):
	import os
	import sys
	import re

	if build == "18":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg18.txt",'r').readlines()
	if build == "19" or build == "1kg":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg19.txt",'r').readlines()

	genes = []
	genes2=[]
	genesnps = []
	genesLong = []
	genesnpsLong = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		if build == "18":
			wsfile = open("/fg/wgas2/rossin/split/SNAPLDFiles/chr"+str(i)+".ws",'r').readlines()
		if build == "19":
			wsfile = open("/fg/wgas2/rossin/split/hg19/chr"+str(i)+".ws",'r').readlines()
		if build == "1kg":
			wsfile = open("/fg/wgas2/rossin/split/1KG/chr"+str(i)+".ws",'r').readlines()					
		for line in wsfile:
			if line.split("\t")[0] in snps:
				count+=1
				snp = line.split("\t")[0]
				chr = "chr"+str(i)
				left = float(line.split("\t")[1])
				right = float(line.strip("\n").split("\t")[2])
				for geneline in genelist:
					geneline = geneline.strip("\n").split("\t")
					if geneline[2]==chr:
						if geneline[3]=="+":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								if not geneline[12] in genes:
									genes.append(geneline[12])
									genesnps.append(line.split("\t")[0])
						if geneline[3]=="-":
							a = float(geneline[4])-regDown
							b = float(geneline[5])+regUp
							if b > left and a < right:
								if not geneline[12] in genes:
									genes.append(geneline[12])
									genesnps.append(line.split("\t")[0])

	ALLgenes = genes[:]
	ALLloci = genesnps[:]

	for i in range(len(genes)):
		if not genes[i] in HC:
			genes[i]="REMOVEME"
			genesnps[i]="REMOVEME"
	PPIgenes = [y for y in genes if y!="REMOVEME"]
	PPIloci = [y for y in genesnps if y!="REMOVEME"]
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	SNPs = list(set(ALLloci))
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	outfile = open(sys.argv[1]+"_MissingGenes",'a')

	for SNP in SNPs:
		outfile.write(str(SNP)+"\t"+str(ALLloci.count(SNP))+"\t"+str(PPIloci.count(SNP))+"\t"+str(round(float(PPIloci.count(SNP))/float(ALLloci.count(SNP)),2))+"\t")
		for i in range(len(PPIgenes)):
			if(PPIloci[i]==SNP):
				outfile.write(str(PPIgenes[i])+",")
		if not SNP in PPIloci:
			outfile.write("-")
		outfile.write("\t")
		find=False
		for i in range(len(ALLgenes)):
			if(ALLloci[i]==SNP and ALLgenes[i] in missedgenes):
				outfile.write(str(ALLgenes[i])+",")
				find=True
		if not find:
			outfile.write("-")
		outfile.write("\n")

	outfile.close()

	results={}
	results[0]=PPIgenes
	results[1]=PPIloci
	return(results)






###Define nearest gene function
def NearestGene(snps,build,HC=HC):
	import os
	import sys
	import re
	if build == "18":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg18.txt",'r').readlines()
	if build == "19" or build == "1kg":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg19.txt",'r').readlines()
	from math import fabs
	outfile = open(sys.argv[1]+"_MissingGenes",'a')
	genes = []
	genesnps = []
	count = 0
	candidateNames=[]
	for i in range(1,23):
		if count==len(snps):
			break
		if build == "18":
			posfile = open("/fg/wgas2/rossin/split/chr"+str(i)+".BIM",'r').readlines()
		if build == "19":
			posfile = open("/fg/wgas2/rossin/split/hg19/chr"+str(i)+".BIM",'r').readlines()
		if build == "1kg":
			posfile = open("/fg/debakker_scratch/ripke/hapmap_ref/impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/eur/ld_info/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1_chr"+str(i)+".impute.phased.bgl.eur.bim",'r').readlines()
		for line in posfile:
			if line.split()[1] in snps:
				distances=[]
				candidates=[]
				candidateNames=[]
				count+=1
				snp = line.split()[1]
				chr = "chr"+str(i)
				pos = float(line.split()[3])
				for geneline in genelist:
					if not geneline.split()[0]=='NA':
						geneline = geneline.strip("\n").split("\t")
						if geneline[2]==chr:
							dist = min([fabs(pos-float(geneline[4])),fabs(pos-float(geneline[5]))])
							distances.append(dist)
							candidates.append(geneline[12])
				candidate = candidates[distances.index(min(distances))]
				if candidate in HC:
					outfile.write(snp+"\t1\t1\t1.0\t"+candidates[distances.index(min(distances))]+"\t-\n")
					genes.append(candidate)
					genesnps.append(line.split()[1])
				elif not candidate in HC:
					outfile.write(snp+"\t1\t0\t0.0\t-\t"+candidates[distances.index(min(distances))]+"\n")

	results={}
	results[0]=genes
	results[1]=genesnps
	return(results)


###Define code to check wether an ensemble ID is a protein
proteinLookup = {}
file = open("/fg/wgas2/rossin/UCSCGeneInfo/EnsemblBiotypes",'r')
file.readline()
for line in file:
	line=line.strip("\n").split(",")
	if len(line)==2:
		proteinLookup[line[0]]=line[1]

#### Region to Gene Function ####
#regions are readlines() from an input file

def RegionToGene(regions,regUp,regDown,build):
	import os
	import sys
	import re

	if build == "18":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg18.txt",'r').readlines()
	elif build == "19" or build == "1kg":
		genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/refseqgenes.hg19.txt",'r').readlines()

	genes = []
	genesnps = []
	genesLong = []
	genesnpsLong = []
	allchrs = []
	allleft = []
	allright = []
	SNPs = []

	for line in regions:
		if not isRegion(line):
			continue
		line = line.rstrip()
		snp = "R"+line.split()[0]
		chr = "chr"+line.split()[1]
		left = float(line.split()[2])
		right = float(line.strip("\n").split()[3])
		for geneline in genelist:
			geneline = geneline.strip("\n").split("\t")
			#if proteinLookup.has_key(geneline[1]):
			#	if not proteinLookup[geneline[1]] == "protein_coding":
			#		continue
			if geneline[2]==chr:
				if geneline[3]=="+":
					a = float(geneline[4])-regUp
					b = float(geneline[5])+regDown
					if b > left and a < right:
						if not geneline[12] in genes:
							genes.append(geneline[12])
							genesnps.append(snp)
							allchrs.append(chr)
							allleft.append(left)
							allright.append(right)
				if geneline[3]=="-":
					a = float(geneline[4])-regDown
					b = float(geneline[5])+regUp
					if b > left and a < right:
						if not geneline[12] in genes:
							genes.append(geneline[12])
							genesnps.append(snp)
							allchrs.append(chr)
							allleft.append(left)
							allright.append(right)
	ALLgenes = genes[:]
	ALLloci = genesnps[:]
	for i in range(len(genes)):
		if not genes[i] in HC:
			genes[i]="REMOVEME"
			genesnps[i]="REMOVEME"
	PPIgenes = [y for y in genes if y!="REMOVEME"]
	PPIloci = [y for y in genesnps if y!="REMOVEME"]
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	outfile = open(sys.argv[1]+"_MissingGenes",'a')
	SNP = ""
	for i in range(len(ALLloci)):
		if ALLloci[i]==SNP:
			continue
		SNP = ALLloci[i]
		outfile.write(str(SNP)+"\t"+str(ALLloci.count(SNP))+"\t"+str(PPIloci.count(SNP))+"\t"+str(float(PPIloci.count(SNP))/float(ALLloci.count(SNP)))+"\t")
		for j in range(len(PPIgenes)):
			if(PPIloci[j]==SNP):
				outfile.write(str(PPIgenes[j])+",")
		outfile.write("\t")
		for j in range(len(ALLgenes)):
			if(ALLloci[j]==SNP and ALLgenes[j] in missedgenes):
				outfile.write(str(ALLgenes[j])+",")
		outfile.write("\n")

	outfile.close()

	results={}
	results[0]=PPIgenes
	results[1]=PPIloci
	return(results)


###Define missing genes function###

###Define region merge function for regions with overlapping genes###

def mergeRegions(seeds,snps):
	import copy
	seedRegion = {}
	regionSeed = {}
	for i in range(len(seeds)):
		if not seedRegion.has_key(seeds[i]):
			seedRegion[seeds[i]]=[]
		seedRegion[seeds[i]].append(snps[i])
	for i in range(len(snps)):
		if not regionSeed.has_key(snps[i]):
			regionSeed[snps[i]]=[]
		regionSeed[snps[i]].append(seeds[i])
	#for each region, figure out overlapping regions
	regionsToJoin = {}
   	for region in regionSeed.keys():
		regionsToJoin[region]=[]
		for otherRegion in regionSeed.keys():
			if otherRegion!=region:
				joint = list(set(regionSeed[region]+regionSeed[otherRegion]))
				if len(joint) < len(regionSeed[region]+regionSeed[otherRegion]):
					regionsToJoin[region].append(otherRegion)
	for region in regionsToJoin.keys():
		for otherRegion in regionsToJoin.keys():
			if region!=otherRegion:
				if region in regionsToJoin[otherRegion]:
					regionsToJoin[region] = regionsToJoin[region]+regionsToJoin[otherRegion]
		regionsToJoin[region] = list(set(regionsToJoin[region]))
	snpsOut=[]
	#for each seed, figure out joint regions
	seeds = seedRegion.keys()
	for seed in seedRegion.keys():
		regions = regionsToJoin[seedRegion[seed][0]]
		regions.sort()
		if len(regions) == 0:
			snpsOut = snpsOut+seedRegion[seed]
		else:
			out = "_".join(regions)
			snpsOut.append(out)
	return([seeds,snpsOut])

def isRegion(input):
	input = input.strip("\n").split()
	if len(input)!=4:
		out=False
	elif input[1].isdigit():
		if len(input)==4 and int(input[1]) in range(1,23) and input[2].isdigit() and input[3].isdigit():
			out=True
		else:
			out=False
	else:
		out=False
	return(out)

def isSnp(input):
	input = input.strip("\n").split()
	if len(input)==1 and (input[0][0:2]=="rs" or input[0][0:3]=="chr"):
		out=True
	else:
		out=False
	return(out)

def isGene(input):
	input = input.strip("\n").split()
	if len(input)==1 and input[0][0:2]!="rs":
		out = True
	else:
		out = False
	return(out)

def isGeneRegion(input):
	input = input.strip("\n").split()
	if len(input)==2:
		out=True
	else:
		out=False
	return(out)


def anyRegion(inputlist):
	for line in inputlist:
		input = line.strip("\n").split()
		if len(input)!=4:
			out=False
			break
		if input[1].isdigit():
			if len(input)==4 and int(input[1]) in range(1,23) and input[2].isdigit() and input[3].isdigit():
				out=True
				break
			else:
				out=False
		else:
			out = False
	return(out)

def anySnp(inputlist):
	for line in inputlist:
		input = line.strip("\n").split()
		if len(input)==1 and (input[0][0:2]=="rs" or input[0][0:3]=="chr"):
			out=True
			break
		else:
			out=False
	return(out)

def anyGene(inputlist):
	for line in inputlist:
		input = line.strip("\n").split()
		if len(input)==1 and input[0][0:2]!="rs" and input[0][0:3]!="chr":
			out=True
			break
		else:
			out=False
	return(out)

def anyGeneRegion(inputlist):
	for line in inputlist:
		input = line.strip("\n").split()
		if len(input)==2:
			out=True
			break
		else:
			out=False
	return(out)
