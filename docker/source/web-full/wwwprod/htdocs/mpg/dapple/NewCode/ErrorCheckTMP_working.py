import os
import sys
import re
import random
import copy

from Functions_slowTMP import *

###Define code to check wether an ensemble ID is a protein
proteinLookup = {}
file = open("/web/dev/htdocs/mpg/dapple/EnsemblBiotypes",'r')
file.readline()
for line in file:
	line=line.strip("\n").split(",")
	if len(line)==2:
		proteinLookup[line[0]]=line[1]

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
build = inputsFile[12].strip("\n")

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
inputlist = list(set(inputlist))

##remove \r \t etc from seedfile
outfile = open("Seeds",'w')
for line in inputlist:
	outfile.write(line.strip()+"\n")
outfile.close()
inputlist = open(seedFile,'r').readlines()

##If they chose color by p-value but did not permute
if plotP=="true" and permute==0:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("You chose to color proteins by DAPPLE p-value but did not choose to permute. You must permute to determine DAPPLE p-values.")
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

##Check that input is the approriate length
if len(inputlist)<2 or len(inputlist)>2500:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Length of input must  be in the range 2-2000.")
	outfile.close()
	sys.exit()	

##Create gene name conversion dictionary
genenames = open("/home/radon00/rossin/DALY/PPI/NewCode/IWtoHugo",'r').readlines()
geneNameConvert={}
for name in genenames:
	name = name.strip("\n").split("\t")
	geneNameConvert[int(name[0])] = name[1]
	geneNameConvert[name[1]]=int(name[0])

inwebProteins = geneNameConvert.keys()

##Figure out what doesn't match in inputlist
print "checking snps..." 
if anySnp(inputlist):
	if build == "18":
		os.system("awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/SNAPLDFILES/*.ws > SNPcheck")
	elif build == "19":
		os.system("awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/hg19/*.ws > SNPcheck")
	elif build == "1kg":
		os.system("awk 'NR==FNR{s[$1]=1;next}{if (s[$1]==1) print }' Seeds /fg/wgas2/rossin/split/1KG/*.ws > SNPcheck")		
	file = open("SNPcheck",'r')
	matchedSNPs = []
	for line in file:
		line = line.split("\t")[0]
		matchedSNPs.append(line)

print "snps grepped."
matched = []
for line in inputlist:
	if isGene(line) or isGeneRegion(line):
		if line.strip("\n").split()[0] in inwebProteins:
			matched.append(inputlist.index(line))
	elif isSnp(line):
		if line.strip("\n") in matchedSNPs:
			matched.append(inputlist.index(line))
	elif isRegion(line):
		matched.append(inputlist.index(line))

for i in range(len(inputlist)):
	if not i in matched:
		warnings.append(inputlist[i].strip("\n"))


if float(len(matched))/float(len(inputlist)) < 0.1:
	outfile = open("ErrorWarningFile.txt",'w')
	os.system("chmod 777 ErrorWarningFile.txt")
	outfile.write("error\n")
	outfile.write("Less than 10% of inputs did not match our system. Please check that format and IDs are correct.")
	outfile.write("The following inputs could not be found:<BR>"+"<BR>".join(warnings))
	outfile.close()
	sys.exit()


###If you made it this far, there are only warnings

outfile = open("ErrorWarningFile.txt",'w')
os.system("chmod 777 ErrorWarningFile.txt")

if len(warnings)>0:
	outfile.write("warning\n")
	outfile.write("The following inputs could not be found:<BR>"+"<BR>".join(warnings))
	outfile.close()
else:
	outfile.write("none\n")
	outfile.close()


