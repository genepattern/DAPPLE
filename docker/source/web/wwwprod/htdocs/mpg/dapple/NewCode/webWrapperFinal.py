import os
import sys

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

##Find input file
seedLMfiles = os.listdir("./")
for file in seedLMfiles:
	if "seedLocusMapping" in file:
		break

##run code
if not genestospecify=='false':
	os.system("python /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py "+file+" permute=0 CIcutoff="+str(CIcutoff)+" plot="+str(plot)+" iterate="+iterate+" compare=true"+" specify="+genestospecify+" keyword="+keyword+" dir="+dir+" nearestgene="+nearestgene+" regulatory="+regulatory+" plotP="+plotP+" collapseCI="+collapseCI+" genome="+genome)
else:
	os.system("python /web/wwwprod/htdocs/mpg/dapple/NewCode/BuildNetwork_slowTMP.py "+file+" permute=0 CIcutoff="+str(CIcutoff)+" plot="+str(plot)+" iterate="+iterate+" compare=true keyword="+keyword+" dir="+dir+" nearestgene="+nearestgene+" regulatory="+regulatory+" plotP="+plotP+" collapseCI="+collapseCI+" genome="+genome)

os.system("chmod 777 ./*")

######email Function#############

SENDMAIL = "/usr/sbin/sendmail" # sendmail location
import os
p = os.popen("%s -t" % SENDMAIL, "w")
p.write("To: "+email+"\n")
p.write("Subject: DAPPLE results\n")
p.write("\n") # blank line separating headers from body
p.write("Your DAPPLE run has completed.\n\n")
p.write("To retrieve your results, return to your job status page (www.broadinstitute.org/mpg/dapple/statusTMP.php?jid="+jobid+") or click on the link below\n\n")
p.write("www.broadinstitute.org/mpg/dapple/results/"+jobid+"/"+keyword+"_results.zip\n\n")
p.write("NOTE: These results are deleted after 7 days\n\n\n")
file = open(keyword+"_directConnections",'r')
p.write("The following direct connections exist between seed proteins from different loci:\n\n")
for line in file:
	p.write(line)
p.write("\n\nThe network parameter values are as follows:\n")
file = open(keyword+"_summary",'r')
for line in file:
	p.write(line)

if permute > 0:
	p.write("\nThe network parameter p-values are as follows:\n")
	file = open(keyword+"_NetStats",'r')
	for line in file:
		p.write(line)
	file = open(keyword+"_seedScores",'r')
	p.write("\nThe seed p-values are as follows:\n")
	for line in file:
		p.write(line)

else:
	p.write("\nYou chose not to permute, so there are no p-values for your network.")

sts = p.close()

os.system("zip "+keyword+"_results.zip "+keyword+"_*"+" Seeds "+keyword+"_GenesToPrioritize")
os.system("chmod 777 "+keyword+"_results.zip")

