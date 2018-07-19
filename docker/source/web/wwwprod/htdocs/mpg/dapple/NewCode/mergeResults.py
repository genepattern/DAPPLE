import os
import sys
import math
import re

inputsFile = open(sys.argv[1],'r').readlines()
seedFile=inputsFile[0].strip("\n")
email=inputsFile[1].strip("\n")
permute=int(inputsFile[2].strip("\n"))
genestospecify=inputsFile[3].strip("\n")
CIcutoff=inputsFile[4].strip("\n")
plot=inputsFile[5].strip("\n")
iterate=inputsFile[6].strip("\n")
keyword=inputsFile[7].strip("\n")
nearestgene = inputsFile[8].strip("\n")
regulatory = inputsFile[9].strip("\n")
plotP = inputsFile[10].strip("\n")
collapseCI = inputsFile[11].strip("\n")
genome = inputsFile[12].strip("\n")

nbreaks=int(sys.argv[2])

##Figure out working directory
tmpdir = seedFile.split("/")
dir = "/"
for component in tmpdir:
	if component == "":
		continue
	if component == tmpdir[len(tmpdir)-1]:
		break
	dir = dir+component+"/"
os.chdir(dir)

##get job ID
jobid = open("statusUpdate",'r').readlines()[0].strip("\n")

##pause while jobs are still running
os.system("bjobs -u dapple -w > tmpbjobs2")
jobs = open("tmpbjobs2",'r').readlines()
stop=False
while not stop:
	os.system("bjobs -u dapple -w > tmpbjobs2")
	jobs = open("tmpbjobs2",'r').readlines()
	stop=True
	for job in jobs:
		for i in range(1,nbreaks):
			if jobid+"_"+str(i) in job:
				stop=False
				break



##concatenate individual files
os.system("cat *_permutedDirectEdgeCount > tmp")
os.system("mv tmp "+keyword+"_permutedDirectEdgeCount")
os.system("cat *_permutedSeedDirectDegreesMean > tmp")
os.system("mv tmp "+keyword+"_permutedSeedDirectDegreesMean")
os.system("cat *_permutedIndirectDegreesMean > tmp")
os.system("mv tmp "+keyword+"_permutedIndirectDegreesMean")
os.system("cat *_permutedCIdegreesMeanPermut > tmp")
os.system("mv tmp "+keyword+"_permutedCIdegreesMeanPermut")


##Write out _NetStats file
import numpy

##Open results file
os.system("cat *_NetStats > all_NetStats")
outfile = open(keyword+"_NetStats",'w')
outfile.write("PARAMETER\tOBSERVED\tEXPECTED\tP_VALUE\n")
diseaseresults = open("all_NetStats",'r').readlines()

#Direct edge count
file = open(keyword+"_permutedDirectEdgeCount",'r')
random = []
for line in file:
	random.append(float(line.strip("\n")))
x = float(diseaseresults[1].split("\t")[1])
outfile.write("Direct Edge Count\t"+str(x)+"\t"+str(numpy.mean(random))+"\t"+str(float(len([i for i in random if i>=x])+1)/float(len(random)+1))+"\n")
print("Direct Edge Count\t"+str(x)+"\t"+str(numpy.mean(random))+"\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

#Direct seed degree mean
file = open(keyword+"_permutedSeedDirectDegreesMean",'r')
random = []
for line in file:
	random.append(float(line.strip("\n")))
x = float(diseaseresults[2].split("\t")[1])
outfile.write("Seed Direct Degrees Mean\t"+str(x)+"\t"+str(numpy.mean(random))+"\t"+str(float(len([i for i in random if i>=x])+1)/float(len(random)+1))+"\n")

#Indirect degrees mean
file = open(keyword+"_permutedIndirectDegreesMean",'r')
random = []
for line in file:
	random.append(float(line.strip("\n")))
x = float(diseaseresults[3].split("\t")[1])
outfile.write("Seed Indirect Degrees Mean\t"+str(x)+"\t"+str(numpy.mean(random))+"\t"+str(float(len([i for i in random if i>=x])+1)/float(len(random)+1))+"\n")

#CI degrees mean
file = open(keyword+"_permutedCIdegreesMeanPermut",'r')
random = []
for line in file:
	random.append(float(line.strip("\n")))
x = float(diseaseresults[4].split("\t")[1])
outfile.write("Seed CI Degrees Mean\t"+str(x)+"\t"+str(numpy.mean(random))+"\t"+str(float(len([i for i in random if i>=x])+1)/float(len(random)+1))+"\n")

print "NetStats merged"

##Write out _SeedScores file
geneSnp = {}
geneCB1 = {}
geneCB2 = {}
successfiles = os.listdir("./")
actualbreaks = 0
for i in range(1,nbreaks):
	if keyword+str(i)+"_seedScores" in successfiles:
		actualbreaks+=1
		file = open(keyword+str(i)+"_seedScores",'r')
		file.readline()
		for line in file:
			line = line.split()
			if line[1]=="NA":
				geneSnp[line[0]]="NA"
				continue
			geneSnp[line[0]]=line[1]
			geneCB1[line[0]]=float(line[4])
			geneCB2[line[0]]=float(line[5])

permute=math.ceil(permute/nbreaks)
outfile = open(keyword+"_seedScores",'w')
outfile.write("GENE\tREGION\tP_uncorrected\tP_corrected\n")
scoreDic = {}
genePermuteCounts = {}
for i in range(1,nbreaks):
	if keyword+str(i)+"_seedPermuteCounts" in successfiles:
		file = open(keyword+str(i)+"_seedPermuteCounts",'r')
		for line in file:
			line = line.strip("\n").split()
			gene = line[0]
			if not genePermuteCounts.has_key(line[0]):
				genePermuteCounts[gene]=[0,0]
			else:
				line[1]=re.sub("NA","0",line[1])
				line[2]=re.sub("NA","0",line[2])
				genePermuteCounts[gene][0]+=int(line[1])
				genePermuteCounts[gene][1]+=int(line[2])

for gene in geneSnp.keys():
	minCount = min([int(genePermuteCounts[gene][0]),int(genePermuteCounts[gene][1])])
	p0=(minCount+1)/(permute*actualbreaks+1)
	if not geneCB1.has_key(gene):
		scoreDic[gene] = [geneSnp[gene],"NA","NA"]
	else:
		p1 = 1-(1-p0)**geneCB1[gene]
		p2 = 1-(1-p1)**geneCB2[gene]
		scoreDic[gene] = [geneSnp[gene],p1,p2]

for key, value in sorted(scoreDic.items(), key=lambda e: e[1][2]):
	outfile.write(key+"\t"+str(value[0])+"\t"+str(value[1])+"\t"+str(value[2])+"\n")
	
outfile.close()

print "SeedScores merged"

##Genes to prioritize
outfile = open(keyword+"_GenesToPrioritize",'w')
for gene in scoreDic.keys():
	if scoreDic[gene][1]!="NA":
		if float(scoreDic[gene][2])<0.05:
			outfile.write(gene+"\n")

outfile.close()

print "Genes prioritized"

##Write CI score file
for file in successfiles:
	if "CIscores" in file:
		CIscores1 = open(file,'r')
		CIscores1.readline()
		correctby = len(CIscores1.readlines())
		CIscores1 = open(file,'r')
		CIscores1.readline()
		scoreDic = {}
		break

outfile = open(keyword+"_CIscores",'w')
outfile.write("PROTEIN\tNUM_BINDERS\tP_VALUE\tP_VALUE_CORRECTED\n")

for line in CIscores1:
	line = line.split()
	gene = line[0]
	os.system("grep -w "+gene+" *_CIscores > tmpCIscores")
	n = float(line[2])*(permute+1)-1
	CIscoresi = open("tmpCIscores",'r')
	for line2 in CIscoresi:
		line2 = line2.split()
		n += float(line2[2])*(permute+1)-1
	p=(n+1)/(permute*actualbreaks+1)
	scoreDic[line[0]]=[line[1],str(p),str(1-(1-p)**correctby)]

for key, value in sorted(scoreDic.items(), key=lambda e: e[1][2]):
	outfile.write(key+"\t"+str(value[0])+"\t"+str(value[1])+"\t"+str(value[2])+"\n")

outfile.close()
