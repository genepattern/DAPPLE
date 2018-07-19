import os
import sys
import math

inputsFile = open(sys.argv[1],'r').readlines()
seedFile=inputsFile[0].strip("\n")
email=inputsFile[1].strip("\n")
permute=int(inputsFile[2].strip("\n"))
genestospecify=inputsFile[3].strip("\n")
CIcutoff=inputsFile[4].strip("\n")
plot=inputsFile[5].strip("\n")
iterate=inputsFile[6].strip("\n")
keyword=inputsFile[7].strip("\n")+sys.argv[2]
nearestgene = inputsFile[8].strip("\n")
regulatory = inputsFile[9].strip("\n")
plotP = inputsFile[10].strip("\n")
collapseCI = inputsFile[11].strip("\n")
genome = inputsFile[12].strip("\n")
zoomedGenes = inputsFile[13].strip("\n")

nbreaks=int(sys.argv[3])
originalKeyword = inputsFile[7].strip("\n")

##Parellelize # permutes
permute = int(math.ceil(permute/nbreaks))
iterate='false'

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

##Write zoomed genes file
out = open("zoomedGeneFile",'w')
out.write("\n".join(zoomedGenes.split(",")))
out.close()

##run code
if not genestospecify=='false':
	os.system("python /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py "+originalKeyword+"1_seedLocusMapping permute="+str(permute)+" CIcutoff="+str(CIcutoff)+" specify="+genestospecify+" keyword="+keyword+" dir="+dir+" nearestgene="+nearestgene+" regulatory="+regulatory+" genome="+genome)
else:
	os.system("python /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py "+originalKeyword+"1_seedLocusMapping permute="+str(permute)+" CIcutoff="+str(CIcutoff)+" keyword="+keyword+" dir="+dir+" nearestgene="+nearestgene+" regulatory="+regulatory+" genome="+genome)
