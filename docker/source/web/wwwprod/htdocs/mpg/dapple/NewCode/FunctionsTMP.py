def buildNet(ppifile,seedsCopy,genesnpCopy,CIcutoff,keyword):
	import copy
	import subprocess
	p=subprocess.Popen(["grep", "-w", "-f", "seedList", ppifile],stdout=subprocess.PIPE)
	ppifile = p.communicate()[0].split("\n")
	ppifile = ppifile[0:len(ppifile)-1]
	seedSeed = {}
	for seed in seedsCopy:
		seedSeed[seed] = []

	seedCI = {}
	for seed in seedsCopy:
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
	seedscopy = copy.deepcopy(seedsCopy)
	nodes=[]
	for line in ppifile:
		line = line.strip("\n").split(d)
		line[0] = int(line[0])
		line[1] = int(line[1])
		nodes.append(line[0])
		nodes.append(line[1])
		if line[0] in seedsCopy and line[1] in seedsCopy:
			seedSeed[line[0]].append(line[1])
			seedCI[line[0]].append(line[1])
		else:
			if line[0] in seedsCopy:
				seedCI[line[0]].append(line[1])
			else:
				if not CIseed.has_key(line[0]):
					CIseed[line[0]] = []
				CIseed[line[0]].append(line[1])
	for seed in seedscopy:
		if not seed in nodes:
			seedsCopy.remove(seed)			
			del genesnpCopy[seed]
##end time consuming part

##Remove seedSeed same locus binding

	seedSeedcopy = copy.deepcopy(seedSeed)
	for seed in seedSeedcopy.keys():
		for otherSeed in seedSeedcopy[seed]:
			if genesnpCopy[seed]==genesnpCopy[otherSeed]:
				seedSeed[seed].remove(otherSeed)

##Count the number of subgraphs in direct network
##Remove any CI that binds only proteins from 1 locus or violates CIcutoff
	CItoDelete = []
	CIseedcopy = copy.deepcopy(CIseed)
	for ci in CIseedcopy.keys():
		cisnp = []
		for seed in CIseed[ci]:
			cisnp.append(genesnpCopy[seed])
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
	for seed in seedsCopy:
		seedDirectDegrees[seed] = len(seedSeed[seed])

	seedIndirectDegrees = {}
	seedSeedIndirect = []
	for seed in seedsCopy:
		numCIconnections = {}
		for ci in seedCI[seed]:
			if CIseed.has_key(ci):
				for otherSeed in CIseed[ci]:
					if otherSeed != seed and genesnpCopy[seed]!=genesnpCopy[otherSeed]:
						if not numCIconnections.has_key(otherSeed):
							numCIconnections[otherSeed] = 0
						numCIconnections[otherSeed] += 1
						indirectPair = [seed,otherSeed]
						indirectPair.sort()
						if not indirectPair in seedSeedIndirect:
							seedSeedIndirect.append(indirectPair)
		numSNPconnections = {}
		for otherSeed in numCIconnections.keys():
			if not numSNPconnections.has_key(genesnpCopy[otherSeed]):
				numSNPconnections[genesnpCopy[otherSeed]] = 0
			if numSNPconnections[genesnpCopy[otherSeed]] < numCIconnections[otherSeed]:
				numSNPconnections[genesnpCopy[otherSeed]] = numCIconnections[otherSeed]
		seedIndirectDegrees[seed] = sum(numSNPconnections.values())

	CIdegreesloc = {}
	CIdegrees = {}
	for ci in CIseed.keys():
		CIdegrees[ci]=len(CIseed[ci])
		cisnp=[]
		for seed in CIseed[ci]:
			if not genesnpCopy[seed] in cisnp:
				cisnp.append(genesnpCopy[seed])
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
								genes.append(int(re.sub(".p","",geneline[0])))
								genesnps.append(line.split("\t")[0])
						if geneline[3]=="-":
							a = float(geneline[4])-regUp
							b = float(geneline[5])+regDown
							if b > left and a < right:
								genes.append(int(re.sub(".p","",geneline[0])))
								genesnps.append(line.split("\t")[0])
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
		HC.append(int(re.sub(".p","",line.split()[0])))
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
					if not geneline.split()[0]=='NA':
						geneline = geneline.strip("\n").split("\t")
						if geneline[2]==chr:
							dist = min([fabs(pos-float(geneline[4])),fabs(pos-float(geneline[5]))])
							distances.append(dist)
							candidates.append(int(re.sub(".p","",geneline[0])))
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

def RegionToGene(regions,regUp=50000,regDown=50000):
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

	for line in regions:
		snp = line.split()[0]
		chr = line.split()[1]
		allchrs.append(chr)
		left = float(line.split()[2])
		allleft.append(left)
		right = float(line.strip("\n").split()[3])
		allright.append(right)
		for geneline in genelist:
			geneline = geneline.strip("\n").split("\t")
			if geneline[2]==chr:
				if geneline[3]=="+":
					a = float(geneline[4])-regUp
					b = float(geneline[5])+regDown
					if b > left and a < right:
						genes.append(int(re.sub(".p","",geneline[0])))
						genesnps.append(line.split()[0])
						genes2.append(geneline[11])
				if geneline[3]=="-":
					a = float(geneline[4])-regDown
					b = float(geneline[5])+regUp
					if b > left and a < right:
						genes.append(int(re.sub(".p","",geneline[0])))
						genesnps.append(line.split()[0])
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
						genesnpsLong.append(line.split("\t")[0])
				if geneline[3]=="-":
					a = float(geneline[4])-regDown
					b = float(geneline[5])+regUp
					if b > left and a < right:
						if not geneline[11] == 'NA':
							genesLong.append(geneline[11])
						else:
							genesLong.append(geneline[1])
						genesnpsLong.append(line.split("\t")[0])
	PPIgenes=genes2
	PPIloci=genesnps

	ALLgenes = genesLong
	ALLloci = genesnpsLong

	SNPs = list(set(ALLloci))
	missedgenes = list(set(ALLgenes)-set(PPIgenes))

	outfile = open(sys.argv[1]+"_MissingGenes",'w')
	outfile.write("REGION\tCHR\tLEFT\tRIGHT\tALLCount\tPPICount\tPercent\tPPIGenes\tMissingGenes\n")
	for i in range(len(SNPs)):
		SNP = SNPs[i]
		outfile.write(str(SNP)+"\t"+str(allchrs[i])+"\t"+str(allleft[i])+"\t"+str(allright[i])+"\t"+str(ALLloci.count(SNP))+"\t"+str(PPIloci.count(SNP))+"\t"+str(float(PPIloci.count(SNP))/float(ALLloci.count(SNP)))+"\t")
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


###Define missing genes function###

def SNPtoGeneCompare(snps,regUp=50000,regDown=50000):
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
								genes.append(int(re.sub(".p","",geneline[0])))
								genesnps.append(line.split("\t")[0])
								genes2.append(geneline[11])
						if geneline[3]=="-":
							a = float(geneline[4])-regDown
							b = float(geneline[5])+regUp
							if b > left and a < right:
								genes.append(int(re.sub(".p","",geneline[0])))
								genesnps.append(line.split("\t")[0])
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
								genesnpsLong.append(line.split("\t")[0])
						if geneline[3]=="-":
							a = float(geneline[4])-regDown
							b = float(geneline[5])+regUp
							if b > left and a < right:
								if not geneline[11] == 'NA':
									genesLong.append(geneline[11])
								else:
									genesLong.append(geneline[1])
								genesnpsLong.append(line.split("\t")[0])

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
