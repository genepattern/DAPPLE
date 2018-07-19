def buildNet(ppifile,dataAnnotFile,seeds,genesnp,CIcutoff,keyword):
	import copy
	ppifile = open(ppifile,'r').readlines()
	ppiannot = open(dataAnnotFile,'r').readlines()
	print "Data loaded."
	dataLoc = {}
	for line in ppiannot:
		line = line.strip("\n").split("\t")
		dataLoc[line[0]] = (int(line[1]),int(line[2]))
	seedFind = {}
	for line in ppiannot:
		line = line.strip("\n").split("\t")
		seedFind[line[0]]=int(line[1])
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
			pair = (string_pair[0], string_pair[1])
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
	print "Dictionaries loaded."
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

	print "Disease parameters calculated."
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


def SNPtoGene(snps,regUp,regDown):
	import os
	import sys
	import re

	genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genes = []
	genesnps = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		wsfile = open("/fg/wgas2/rossin/split/SNAPLDFiles/chr"+str(i)+".ws",'r').readlines()
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
								genes.append(geneline[11])
								genesnps.append(line.split("\t")[11])
						if geneline[3]=="-":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								genes.append(geneline[11])
								genesnps.append(line.split("\t")[11])
	results={}
	results[0]=genes
	results[1]=genesnps
	return(results)

###Define nearest gene function

def NearestGene(snps):
	import os
	import sys
	import re
	HC = []
	genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genelist.remove('IW\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\thugo\tgene.score\trank.score\tind.SNPs\tcorrected\trank.corrected\n')
	for line in genelist:
		HC.append([11])
	from math import fabs
	genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ens_scored.txt",'r').readlines()
	outfile = open(sys.argv[1]+"_MissingGenes",'w')
	outfile.write("SNP\tCLOSESTGENEPPI\tCLOSEGENEMISSED\n")
	genes = []
	genesnps = []
	count = 0
	candidateNames=[]
	for i in range(1,23):
		if count==len(snps):
			break
		posfile = open("/fg/wgas2/rossin/split/chr"+str(i)+".BIM",'r').readlines()
		for line in posfile:
			if line.split()[1] in snps:
				distances=[]
				candidates=[]
				candidateNames=[]
				count+=1
				snp = line.split()[1]
				chr = str(i)
				pos = float(line.split()[3])
				for geneline in genelist:
					if not geneline.split()[11]=='NA':
						geneline = geneline.strip("\n").split("\t")
						if geneline[2]==chr:
							dist = min([fabs(pos-float(geneline[4])),fabs(pos-float(geneline[5]))])
							distances.append(dist)
							candidates.append(geneline[11])
							candidateNames.append(geneline[11])
				candidate = candidates[distances.index(min(distances))]
				if candidate in HC:
					outfile.write(snp+"\t"+candidateNames[distances.index(min(distances))]+"\t\n")
				elif not candidate in HC:
					outfile.write(snp+"\t\t"+candidateNames[distances.index(min(distances))]+"\n")
				genes.append(candidate)
				genesnps.append(line.split()[1])
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

def RegionToGene(regions,regUp,regDown):
	import os
	import sys
	import re

	genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genelistfull = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ens_scored.txt",'r').readlines()
	genes = []
	genesnps = []
	genes2=[]
	genesLong = []
	genesnpsLong = []
	allchrs = []
	allleft = []
	allright = []
	SNPs = []

	for line in regions:
		line = line.rstrip()
		snp = line.split()[0]
		chr = line.split()[1]
		left = float(line.split()[2])
		right = float(line.strip("\n").split()[3])
		for geneline in genelist:
			geneline = geneline.strip("\n").split("\t")
			if geneline[2]==chr:
				if geneline[3]=="+":
					a = float(geneline[4])-regUp
					b = float(geneline[5])+regDown
					if b > left and a < right:
						genes.append(geneline[11])
						genesnps.append(line.split()[11])
						genes2.append(geneline[11])
				if geneline[3]=="-":
					a = float(geneline[4])-regDown
					b = float(geneline[5])+regUp
					if b > left and a < right:
						genes.append(geneline[11])
						genesnps.append(line.split()[11])
						genes2.append(geneline[11])						
		for geneline in genelistfull:
			geneline = geneline.strip("\n").split("\t")
			if proteinLookup.has_key(geneline[1]):
				if not proteinLookup[geneline[1]] == "protein_coding":
					continue
			if geneline[2]==chr:
				if geneline[3]=="+":
					a = float(geneline[4])-regUp
					b = float(geneline[5])+regDown
					if b > left and a < right:
						if not geneline[11] == 'NA':
							genesLong.append(geneline[11])
						else:
							genesLong.append(geneline[1])
						genesnpsLong.append(line.split()[11])
						SNPs.append(snp)
						allchrs.append(chr)
						allleft.append(left)
						allright.append(right)
				if geneline[3]=="-":
					a = float(geneline[4])-regDown
					b = float(geneline[5])+regUp
					if b > left and a < right:
						if not geneline[11] == 'NA':
							genesLong.append(geneline[11])
						else:
							genesLong.append(geneline[1])
						genesnpsLong.append(line.split()[11])
						SNPs.append(snp)
						allchrs.append(chr)
						allleft.append(left)
						allright.append(right)
	PPIgenes=genes2
	PPIloci=genesnps
	ALLgenes = genesLong
	ALLloci = genesnpsLong
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	outfile = open(sys.argv[1]+"_MissingGenes",'w')
	outfile.write("REGION\tCHR\tLEFT\tRIGHT\tALLCount\tPPICount\tPercent\tPPIGenes\tMissingGenes\n")
	SNP = ""
	for i in range(len(SNPs)):
		if SNPs[i]==SNP:
			continue
		SNP = SNPs[i]
		outfile.write(str(SNP)+"\t"+str(allchrs[i])+"\t"+str(allleft[i])+"\t"+str(allright[i])+"\t"+str(ALLloci.count(SNP))+"\t"+str(PPIloci.count(SNP))+"\t"+str(float(PPIloci.count(SNP))/float(ALLloci.count(SNP)))+"\t")
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
	results[0]=genes
	results[1]=genesnps
	return(results)


###Define missing genes function###

def SNPtoGeneCompare(snps,regUp,regDown):
	import os
	import sys
	import re

	genelist = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genelistfull = open("/fg/wgas2/rossin/UCSCGeneInfo/gene.info_hg18_ens_scored.txt",'r').readlines()
	genes = []
	genes2=[]
	genesnps = []
	genesLong = []
	genesnpsLong = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		wsfile = open("/fg/wgas2/rossin/split/SNAPLDFiles/chr"+str(i)+".ws",'r').readlines()
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
								genes.append(geneline[11])
								genesnps.append(line.split("\t")[11])
								genes2.append(geneline[11])
						if geneline[3]=="-":
							a = float(geneline[4])-regDown
							b = float(geneline[5])+regUp
							if b > left and a < right:
								genes.append(geneline[11])
								genesnps.append(line.split("\t")[11])
								genes2.append(geneline[11])
				for geneline in genelistfull:
					geneline = geneline.strip("\n").split("\t")
					if proteinLookup.has_key(geneline[1]):
						if not proteinLookup[geneline[1]] == "protein_coding":
							continue
					if geneline[2]==chr:
						if geneline[3]=="+":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								if not geneline[11] == 'NA':
									genesLong.append(geneline[11])
								else:
									genesLong.append(geneline[1])
								genesnpsLong.append(line.split("\t")[11])
						if geneline[3]=="-":
							a = float(geneline[4])-regDown
							b = float(geneline[5])+regUp
							if b > left and a < right:
								if not geneline[11] == 'NA':
									genesLong.append(geneline[11])
								else:
									genesLong.append(geneline[1])
								genesnpsLong.append(line.split("\t")[11])

	PPIgenes=genes2
	PPIloci=genesnps

	ALLgenes = genesLong
	ALLloci = genesnpsLong

	SNPs = list(set(ALLloci))
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	outfile = open(sys.argv[1]+"_MissingGenes",'w')
	outfile.write("SNP\tALLCount\tPPICount\tPercent\tPPIGenes\tMissingGenes\n")
	for SNP in SNPs:
		outfile.write(str(SNP)+"\t"+str(ALLloci.count(SNP))+"\t"+str(PPIloci.count(SNP))+"\t"+str(float(PPIloci.count(SNP))/float(ALLloci.count(SNP)))+"\t")
		for i in range(len(PPIgenes)):
			if(PPIloci[i]==SNP):
				outfile.write(str(PPIgenes[i])+",")
		outfile.write("\t")
		for i in range(len(ALLgenes)):
			if(ALLloci[i]==SNP and ALLgenes[i] in missedgenes):
				outfile.write(str(ALLgenes[i])+",")
		outfile.write("\n")

	outfile.close()

	results={}
	results[0]=genes
	results[1]=genesnps
	return(results)




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
			out = ""
			for region in regions:
				out = out + "_" + region
			snpsOut.append(out)
	return([seeds,snpsOut])

def isRegion(input):
	input = input.strip("\n").split()
	if len(input)==4 and int(input[1]) in range(1,23) and input[2].isdigit() and input[3].isdigit():
		return('true')
	else:
		return('false')

def isSNP(input):
	input = input.strip("\n").split()
	if len(input)==1 and input[0][0:2]=="rs":
		return('true')
	else:
		return('false')

def isGene(input):
	input = input.strip("\n").split()
	if len(input)==1 and input[0][0:2]!="rs":
		return('true')
	else:
		return('false')

def isGeneRegion(input):
	input = input.strip("\n").split()
	if len(input)==2:
		return('true')
	else:
		return('false')
