import os
import sys
import re
import random
import copy

#####SNP to Gene Function#######

def SNPtoGene(snps,regUp,regDown):
	import os
	import sys
	import re

	genelist = open("/web/dev/htdocs/mpg/dapple/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genes = []
	genesnps = []

	count = 0
	for i in range(1,23):
		if count==len(snps):
			break
		wsfile = open("/web/dev/htdocs/mpg/dapple/wingspanFiles/chr"+str(i)+".ws",'r').readlines()
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


###Define code to check wether an ensemble ID is a protein
proteinLookup = {}
file = open("/web/dev/htdocs/mpg/dapple/EnsemblBiotypes",'r')
file.readline()
for line in file:
	line=line.strip("\n").split(",")
	if len(line)==2:
		proteinLookup[line[0]]=line[1]

#### Region to Gene Function ####
def RegionToGene(regions,regUp=50000,regDown=50000):
	import os
	import sys
	import re

	genelist = open("/web/dev/htdocs/mpg/dapple/gene.info_hg18_ensHC_scored.txt",'r').readlines()
	genelistfull = open("/web/dev/htdocs/mpg/dapple/gene.info_hg18_ens_scored.txt",'r').readlines()
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


##input can be SNPs, genes or regions
warnings=[]
errors=[]
inputlist = open(seedFile,'r').readlines()

##Check for lines of whitespace
for line in inputlist:
	line = line.rstrip()
	if line=="":
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("Please remove any lines of whitespace from your input.")
		outfile.close()
		sys.exit()

##Check that if they selected nearest gene, that they input SNPs
if nearestgene=="true":
	if not inputlist[0][0:2]=="rs":
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("You must input SNPs if you chose to use the nearest gene.")
		outfile.close()
		sys.exit()
		
##Check that keyword is valid
if "_" in keyword or "-" in keyword or "/" in keyword or "\\" in keyword or "," in keyword:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Invalid description. Please do not use any of the forbidden characters.")
	outfile.close()
	sys.exit()

##Check that the keyword is 1 word
if len(keyword.split()) > 1:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Please make keyword 1 word.")
	outfile.close()
	sys.exit()

##Create gene name conversion dictionary
genenames = open("/home/radon00/rossin/DALY/PPI/NewCode/IWtoHugo",'r').readlines()
geneNameConvert={}
for name in genenames:
	name = name.strip("\n").split("\t")
	geneNameConvert[int(name[0])] = name[1]
	geneNameConvert[name[1]]=int(name[0])

##if input is a list of SNPs
if inputlist[0][0:2]=="rs" and len(inputlist[0].split())!=4:
	###check that there is one snp per line
	if len(inputlist[0].split())>1:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("If inputting SNPs, please only have 1 SNP per line.")
		outfile.close()
		sys.exit()

	###check that there are no more than 150 snps
	if len(inputlist) > 10000:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("The maximum number of SNPs to input is 10000.")
		outfile.close()
		sys.exit()

    ###check that there is more than 1 snp
	if len(inputlist)==1:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("You must input at least 2 snps.")
		outfile.close()
		sys.exit()
		
	###check that all snps begin with "rs":
	for snp in inputlist:
		if snp[0:2]!="rs":
			outfile = open("ErrorWarningFile.txt",'w')
			os.system("chmod 777 ErrorWarningFile.txt")
			outfile.write("error\n")
			outfile.write("All SNPs must begin with 'rs'. Do not mix SNPs and genes.")
			outfile.close()
			sys.exit()
	
	os.system("rm SNPcheck")
	os.system("sed s/\r//g "+seedFile+" > tmp")
	os.system("mv tmp "+seedFile)
	os.system("grep -w -f "+seedFile+" /web/dev/htdocs/mpg/dapple/wingspanFiles/*.ws > SNPcheck")
	
	file = open("SNPcheck",'r')
	matchedSNPs = []
	for line in file:
		line = line.split("\t")[0]
		line = line.split(":")[1]
		matchedSNPs.append(line)		
	snps = []
	for line in inputlist:
		snps.append(line.strip("\n").strip("\r"))
	unmatchedSNPs = list(set(snps)-set(matchedSNPs))
	unmatchedWarning=""
	for SNP in unmatchedSNPs:
		unmatchedWarning = unmatchedWarning+SNP+","
	unmatchedWarning=unmatchedWarning.strip(",")
	if len(unmatchedWarning)>0:
		warnings.append("The following SNPs were not matched: "+unmatchedWarning+"; Consider using proxies")

	x = SNPtoGene(snps,50000,50000)
	seeds = x[0]
	snps = x[1]

##if input is a list of genes
elif (len(inputlist[0].split())==1 and inputlist[0][0:2]!="rs") or len(inputlist[0].split())==2:
	if len(inputlist) > 5000:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("Please input fewer than 5000 entries.")
		outfile.close()
		sys.exit()
	if len(inputlist)==1:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("You must input at least 2 seed genes.")
		outfile.close()
		sys.exit()

	print("Input is a list of genes. Extracting genes...")
	count=0
	seeds=[]
	snps=[]
	unmatchedGenes=""
	for line in inputlist:
		if geneNameConvert.has_key(line.strip("\n").split()[0]):
			seeds.append(geneNameConvert[line.strip("\n").split()[0]])
			if len(inputlist[0].split())==2:
				snps.append(line.strip("\n").split()[1])
			else:
				snps.append(str(count))
				count+=1
		else:
			unmatchedGenes = unmatchedGenes+line.strip("\n").split()[0]+","
	if len(unmatchedGenes)>0:
		warnings.append("The following input genes could not be found: "+unmatchedGenes+"; Make sure you are only entering protein-coding genes and that you are using the correct gene IDs.")

##if input is a list of regions
elif len(inputlist[0].split())==4:
	for line in inputlist:
		line = line.strip("\n").split()
		if not line[1] in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23']:
			outfile = open("ErrorWarningFile.txt",'w')
			os.system("chmod 777 ErrorWarningFile.txt")
			outfile.write("error\n")
			outfile.write("When entering regions, the second column must be a number in [1,23], where 23 is the X chromosome")
			outfile.close()
			sys.exit()
		if not line[2].isdigit() or not line[2].isdigit():
			outfile = open("ErrorWarningFile.txt",'w')
			os.system("chmod 777 ErrorWarningFile.txt")
			outfile.write("error\n")
			outfile.write("When entering regions, the third and fourth columns must be integers")
			outfile.close()
			sys.exit()
	if len(inputlist) > 2000:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("Please do not enter more than 2000 regions.")
		outfile.close()
		sys.exit()
	if len(inputlist) == 1:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("You must input at least 2 regions.")
		outfile.close()
		sys.exit()
	if not inputlist[0].split()[1] in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23']:
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("When entering regions, the second column must be a number in [1,23]")
		outfile.close()
		sys.exit()
	elif not inputlist[0].split()[2].isdigit() or not inputlist[0].split()[3].isdigit():
		outfile = open("ErrorWarningFile.txt",'w')
		os.system("chmod 777 ErrorWarningFile.txt")
		outfile.write("error\n")
		outfile.write("The region boundaries must be integers.")
		outfile.close()
		sys.exit()
	else:
		print("Input is a list of regions. Converting regions to genes...")
		x = RegionToGene(inputlist,50000,50000)
		seeds = x[0]
		snps = x[1]
else:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Input format is invalid. Input must be one of the following:\n-SNP list (one SNP per line, each identified by an rs number)\n-Gene list (each line consisting of one gene or one gene and a region ID, separated by a tab)\n-One region per line (each line consists of a region ID, chromosome, left bound and right bound, all tab separated)")
	outfile.close()
	sys.exit()

##Exit if too few genes match PPI data

if float(len(seeds))/float(len(inputlist)) < 0.1:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Less than 10% of seed genes are in the PPI database and DAPPLE did not run. Check to make sure IDs are Hugo IDs - for example, \"ATXN1\". If your input was a list of SNPs, make sure they are in HapMap release 22.")
	outfile.close()
	sys.exit()

###If you made it this far, there are only warnings

outfile = open("ErrorWarningFile.txt",'w')
os.system("chmod 777 ErrorWarningFile.txt")
if len(warnings)>0:
	outfile.write("warning\n")
	for warning in warnings:
		outfile.write(warning+"\n")
	outfile.close()
else:
	outfile.write("none\n")
	outfile.close()

