import gzip

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
	seedSeed = {}
	for seed in seeds:
		seedSeed[seed] = []

	seedCI = {}
	for seed in seeds:
		seedCI[seed] = []

	CIseed = {}

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
		subData = ppifile[dataLoc[seed][0]:dataLoc[seed][1]]
		for line in subData:
			string_pair = line.strip("\n").split(d)
			pair = (int(string_pair[0]), int(string_pair[1]))
			if pair[1] in seeds:
				seedSeed[pair[0]].append(pair[1])
				seedCI[pair[0]].append(pair[1])
				if not CIseed.has_key(pair[0]):
					CIseed[pair[0]] = []
				CIseed[pair[0]].append(pair[1])
			else:
				seedCI[pair[0]].append(pair[1])
				if not CIseed.has_key(pair[1]):
					CIseed[pair[1]] = []
				CIseed[pair[1]].append(pair[0])

##end time consuming part

##Remove seedSeed same locus binding

	seedSeedcopy = copy.deepcopy(seedSeed)
	for seed in seedSeedcopy.keys():
		for otherSeed in seedSeedcopy[seed]:
			if genesnp[seed]==genesnp[otherSeed]:
				seedSeed[seed].remove(otherSeed)

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

##Remove any CI that binds only proteins from 1 locus or violates CIcutoff

	CItoDelete = []
	CIseedcopy = copy.deepcopy(CIseed)
	for ci in CIseedcopy.keys():
		cisnp = []
		for seed in CIseed[ci]:
			cisnp.append(genesnp[seed])
		if len(list(set(cisnp))) == 1 and CIcutoff > 1:
			del CIseed[ci]
			CItoDelete.append(ci)
		elif len(CIseed[ci]) < CIcutoff:
			del CIseed[ci]
			CItoDelete.append(ci)

	for ci in CItoDelete:
		for seed in seedCI.keys():
			if ci in seedCI[seed]:
				seedCI[seed].remove(ci)

##Create degrees dictionaries

	seedDirectDegrees = {}
	for seed in seeds:
		seedDirectDegrees[seed] = len(seedSeed[seed])

	seedIndirectDegrees = {}
	seedSeedIndirect = []
	for seed in seeds:
		numCIconnections = {}
		for ci in seedCI[seed]:
			if CIseed.has_key(ci):
				for otherSeed in CIseed[ci]:
					if otherSeed != seed and genesnp[seed]!=genesnp[otherSeed]:
						if not numCIconnections.has_key(otherSeed):
							numCIconnections[otherSeed] = 0
						numCIconnections[otherSeed] += 1
						indirectPair = [seed,otherSeed]
						indirectPair.sort()
						if not indirectPair in seedSeedIndirect:
							seedSeedIndirect.append(indirectPair)
		numSNPconnections = {}
		for otherSeed in numCIconnections.keys():
			if not numSNPconnections.has_key(genesnp[otherSeed]):
				numSNPconnections[genesnp[otherSeed]] = 0
			if numSNPconnections[genesnp[otherSeed]] < numCIconnections[otherSeed]:
				numSNPconnections[genesnp[otherSeed]] = numCIconnections[otherSeed]
		seedIndirectDegrees[seed] = sum(numSNPconnections.values())

	CIdegreesloc = {}
	CIdegrees = {}
	for ci in CIseed.keys():
		CIdegrees[ci]=len(CIseed[ci])
		cisnp=[]
		for seed in CIseed[ci]:
			if not genesnp[seed] in cisnp:
				cisnp.append(genesnp[seed])
			CIdegreesloc[ci]=len(cisnp)

##summarize results
	def mean(x):
		while 0 in x:
			x.remove(0)
		if len(x)==0:
			return(0)
		return(float(sum(x))/float(len(x)))

##return results
	results={}
	results["directEdgesCount"] = sum(seedDirectDegrees.values())/2
	results["seedDirectDegrees"] = seedDirectDegrees
	results["seedDirectDegreesMean"] = mean(seedDirectDegrees.values())
	results["seedIndirectDegrees"] = seedIndirectDegrees
	results["seedIndirectDegreesMean"] = mean(seedIndirectDegrees.values())
#	results["CIdegrees"] = CIdegreesloc
#	results["CIdegreesMean"] = mean(CIdegreesloc.values())
	results["CIdegrees"] = CIdegrees
	results["CIdegreesMean"] = mean(CIdegrees.values())
	results["seedSeed"] = seedSeed
	results["seedSeedIndirect"] = seedSeedIndirect
	results["seedCI"] = seedCI
	results["CIseed"] = CIseed
#	results["subGraphs"] = subGraphs
	return(results)

#####SNP to Gene Function#######
###Define list of proteins in the InWeb3 HC set
HC=[]
HCfile = open("/data/InWeb3_HC_Red_hugo",'r').readlines()
for line in HCfile:
	HC.append(line.strip("\n").split()[0])
	HC.append(line.strip("\n").split()[1])
HC = list(set(HC))

def SNPtoGeneOLD(snps,regUp,regDown,build):
	import os
	import sys
	import re
	if build == "18":
		genelist = open("/data/refseq/refseqgenes.hg18.txt",'r').readlines()

	if build == "19" or build == "1kg":
		genelist = open("/data/refseq/refseqgenes.hg19.txt",'r').readlines()
	genes = []
	genesnps = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		if build == "18":
			wsfile = gzip.open("/data/wingspan/hg18/chr"+str(i)+".ws.gz",'r').readlines()
		if build == "19" or build == "1kg":
			wsfile = gzip.open("/data/wingspan/hg19/chr"+str(i)+".ws.gz",'r').readlines()
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
		genelist = open("/data/refseq/refseqgenes.hg18.txt",'r').readlines()
	if build == "19" or build == "1kg":
		genelist = open("/data/refseq/refseqgenes.hg19.txt",'r').readlines()

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
			wsfile = gzip.open("/data/wingspan/hg18/chr"+str(i)+".ws.gz",'r').readlines()
		if build == "19":
			wsfile = gzip.open("/data/wingspan/hg19/chr"+str(i)+".ws.gz",'r').readlines()
		if build == "1kg":
			wsfile = gzip.open("/data/wingspan/1KG/chr"+str(i)+".ws.gz",'r').readlines()

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
		genelist = open("/data/refseq/refseqgenes.hg18.txt",'r').readlines()
	if build == "19" or build == "1kg":
		genelist = open("/data/refseq/refseqgenes.hg19.txt",'r').readlines()
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
			posfile = gzip.open("/data/bim/hg18/chr"+str(i)+".BIM.gz",'r').readlines()
		if build == "19":
			posfile = gzip.open("/data/bim/hg19/chr"+str(i)+".BIM.gz",'r').readlines()
		if build == "1kg":
			#Hope those are the real files, copied from a likely backup, might not be imputed/phased version?
			#posfile = open("/fg/debakker_scratch/ripke/hapmap_ref/impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/eur/ld_info/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1_chr"+str(i)+".impute.phased.bgl.eur.bim",'r').readlines()
			posfile = gzip.open("data/bim/1KG/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1_chr"+str(i)+".eur.bfile.bim.gz",'r').readlines()
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
file = open("/data/EnsemblBiotypes",'r')
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
		genelist = open("/data/refseq/refseqgenes.hg18.txt",'r').readlines()
	elif build == "19" or build == "1kg":
		genelist = open("/data/refseq/refseqgenes.hg19.txt",'r').readlines()

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

def buildNetRandom(ppifile,dataAnnotFile,seeds,genesnp,CIcutoff,keyword,randomSeed):
	import copy
	ppifile = open(ppifile,'r').readlines()
	ppiannot = open(dataAnnotFile,'r').readlines()

	#mix seed and genesnp
	shuffledNodes = getNewMapping(ppifile, randomSeed) # -> line 767
#	shuffledNodesReverse = {y:x for x,y in shuffledNodes.iteritems()}

	seedsDict = {}
	for index, item in enumerate(seeds):
		# print index, item, " -> ", shuffledNodes[item] 
		seeds[index] = shuffledNodes[item]
		seedsDict[shuffledNodes[item]] = 1
		genesnp[shuffledNodes[item]] = genesnp[item]
	# print

	# print genesnp
	for key in genesnp.keys():
		if key not in seedsDict:
			del genesnp[key]
	# print genesnp
	# print

	dataLoc = {}
	for line in ppiannot:
		line = line.strip("\n").split("\t")
		dataLoc[int(line[0])] = (int(line[1]),int(line[2]))
	seedFind = {}

	seedSeed = {}
	for seed in seeds:
		seedSeed[seed] = []

	seedCI = {}
	for seed in seeds:
		seedCI[seed] = []

	CIseed = {}

	if "\t" in ppifile[1]:
		d = "\t"
	elif "-" in ppifile[1]:
		d = "-"
	else:
		print(ppifile[1])
		print(" does not have a delimiter")

##begin time consuming part (rate-limiting)

	seedscopy = copy.deepcopy(seeds)							# afraid of removing keys while iterating
	for seed in seedscopy:										#
		if not dataLoc.has_key(seed):							# remove seeds without representation in network
			seeds.remove(seed)									#
			del genesnp[seed]									#
			continue											#
		subData = ppifile[dataLoc[seed][0]:dataLoc[seed][1]]	# get chunk with interactions from index
		for line in subData:									#
			string_pair = line.strip("\n").split(d)				#
			pair = (int(string_pair[0]), int(string_pair[1]))	# Only shuffle interactors 
			if pair[1] in seeds:								# if interactor in seeds
				seedSeed[pair[0]].append(pair[1])				#   add to seedSeed		seed -> int
				seedCI[pair[0]].append(pair[1])					#   add to seedCI		seed -> int
				if not CIseed.has_key(pair[0]):					#
					CIseed[pair[0]] = []						#
				CIseed[pair[0]].append(pair[1])					#   add to CIseed
			else:												# else
				seedCI[pair[0]].append(pair[1])					#   add to seedCI		seed->int
				if not CIseed.has_key(pair[1]):					#
					CIseed[pair[1]] = []						#
				CIseed[pair[1]].append(pair[0])					#   add to CIseed		int -> seed

##end time consuming part

##Remove seedSeed same locus binding
#	print genesnp
	seedSeedcopy = copy.deepcopy(seedSeed)						# afraid of removing keys while iterating
	for seed in seedSeedcopy.keys():							# for each seed
		for otherSeed in seedSeedcopy[seed]:					#   for each interacting seed
			if genesnp[seed]==genesnp[otherSeed]:				#     if same genesnp, reverse lookup, NOT forward
				seedSeed[seed].remove(otherSeed)				#       remove interactor from genesnp

##Remove any CI that binds only proteins from 1 locus or violates CIcutoff

	CItoDelete = []
	CIseedcopy = copy.deepcopy(CIseed)
	for ci in CIseedcopy.keys():
		cisnp = []
		for seed in CIseed[ci]:
			cisnp.append(genesnp[seed])
		if len(list(set(cisnp))) == 1 and CIcutoff > 1:
			del CIseed[ci]
			CItoDelete.append(ci)
		elif len(CIseed[ci]) < CIcutoff:
			del CIseed[ci]
			CItoDelete.append(ci)

	for ci in CItoDelete:
		for seed in seedCI.keys():
			if ci in seedCI[seed]:
				seedCI[seed].remove(ci)

##Create degrees dictionaries

	seedDirectDegrees = {}
	for seed in seeds:
		seedDirectDegrees[seed] = len(seedSeed[seed])

	seedIndirectDegrees = {}
	seedSeedIndirect = []
	for seed in seeds:
		numCIconnections = {}
		for ci in seedCI[seed]:
			if CIseed.has_key(ci):
				for otherSeed in CIseed[ci]:
					if otherSeed != seed and genesnp[seed]!=genesnp[otherSeed]:
						if not numCIconnections.has_key(otherSeed):
							numCIconnections[otherSeed] = 0
						numCIconnections[otherSeed] += 1
						indirectPair = [seed,otherSeed]
						indirectPair.sort()
						if not indirectPair in seedSeedIndirect:
							seedSeedIndirect.append(indirectPair)
		numSNPconnections = {}
		for otherSeed in numCIconnections.keys():
			if not numSNPconnections.has_key(genesnp[otherSeed]):
				numSNPconnections[genesnp[otherSeed]] = 0
			if numSNPconnections[genesnp[otherSeed]] < numCIconnections[otherSeed]:
				numSNPconnections[genesnp[otherSeed]] = numCIconnections[otherSeed]
		seedIndirectDegrees[seed] = sum(numSNPconnections.values())

	CIdegreesloc = {}
	CIdegrees = {}
	for ci in CIseed.keys():
		CIdegrees[ci]=len(CIseed[ci])
		cisnp=[]
		for seed in CIseed[ci]:
			if not genesnp[seed] in cisnp:
				cisnp.append(genesnp[seed])
			CIdegreesloc[ci]=len(cisnp)

##summarize results
	def mean(x):
		while 0 in x:
			x.remove(0)
		if len(x)==0:
			return(0)
		return(float(sum(x))/float(len(x)))

##return results
	results={}
	results["directEdgesCount"] = sum(seedDirectDegrees.values())/2
	results["seedDirectDegrees"] = seedDirectDegrees
	results["seedDirectDegreesMean"] = mean(seedDirectDegrees.values())
	results["seedIndirectDegrees"] = seedIndirectDegrees
	results["seedIndirectDegreesMean"] = mean(seedIndirectDegrees.values())
	results["CIdegrees"] = CIdegrees
	results["CIdegreesMean"] = mean(CIdegrees.values())
	results["seedSeed"] = seedSeed
	results["seedSeedIndirect"] = seedSeedIndirect
	results["seedCI"] = seedCI
	results["CIseed"] = CIseed
	return(results)

def getNewMapping(ppifile, randomSeed):
	import random
	import sys

#	print >> sys.stderr, ppifileName

#	data = open(ppifileName,'r')
	degree = {}
	pairs = []
	for line in ppifile:
		line = line.strip().split("\t")
		pairs.append([line[0],line[1]])
		if line[0] not in degree:
			degree[line[0]] = 1
		else:
			degree[line[0]]+=1
		if line[1] not in degree:
			degree[line[1]] = 1
		else:
			degree[line[1]]+=1

	degreeOptions = list(set((degree.values())))
	degreeOptions.sort()

	newnodes = {}                                                                           # empty newnodes
	i=0
	while i <= len(degreeOptions)-1:                                                        # while i smaller than max degree
		tmpnodes = []                                                                       #   empty new node labels
		binsize=0                                                                           #   initialize binsize
		while binsize < 20:                                                                 #   while binsize less than 20
			tmpdegree = degreeOptions[i]                                                    #     get current degree
			tmpnodes = tmpnodes+[node for node, dg in degree.items() if dg == tmpdegree]    #     get all nodes with same degree
			binsize = len(tmpnodes)                                                         #     binsize is # of tmp nodes
			i+=1                                                                            #     next degree
			if i >= len(degreeOptions)-1:                                                   #     stop if we are out og genes
				break
		if not randomSeed == "NA":
			random.seed(randomSeed)
		tmpshuffle = random.sample(tmpnodes,len(tmpnodes))                                  #   shuffle node labels
		for j in range(len(tmpnodes)):                                                      #   for all genes in this bin
			newnodes[int(tmpnodes[j])] = int(tmpshuffle[j])                                           #     assign new labels
	return newnodes
