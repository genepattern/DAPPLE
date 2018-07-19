import os
import sys

outfile = open("Heiko_700_NetStats_merge",'w')

numpermute = len(open("heiko_input_700_P62_permutedDirectEdgeCount",'r').readlines())
diseaseresults = open("heiko_input_700_P62_NetStats",'r').readlines()
x = float(diseaseresults[1].split("\t")[1])
os.system("rm heiko_input_700_permutedDirectEdgeCount")
os.system("cat heiko_input_700_P*_permutedDirectEdgeCount > heiko_input_700_permutedDirectEdgeCount")
random=[]
file = open("heiko_input_700_permutedDirectEdgeCount",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Direct Edge Count\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[2].split("\t")[2])
os.system("rm heiko_input_700_permutedSeedDirectDegreesMean")
os.system("cat heiko_input_700_P*_permutedSeedDirectDegreesMean > heiko_input_700_permutedSeedDirectDegreesMean")
random=[]
file = open("heiko_input_700_permutedSeedDirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Direct Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[3].split("\t")[2])
os.system("rm heiko_input_700_permutedIndirectDegreesMean")
os.system("cat heiko_input_700_P*_permutedIndirectDegreesMean > heiko_input_700_permutedIndirectDegreesMean")
random=[]
file = open("heiko_input_700_permutedIndirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Indirect Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[4].split("\t")[2])
os.system("rm heiko_input_700_permutedCIdegreesMeanPermut")
os.system("cat heiko_input_700_P*_permutedCIdegreesMeanPermut > heiko_input_700_permutedCIdegreesMeanPermut")
random=[]
file = open("heiko_input_700_permutedCIdegreesMeanPermut",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("CI Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

outfile.close()

os.system("cat heiko_input_700_P*_SeedScores > heiko_input_700_SeedScores")
file = open("heiko_input_700_SeedScores",'r')
genes = {}
for line in file:
	line = line.strip("\n").split()
	if line[0]=="GENE" or line[2]=="NA":
		continue
	if not line[0] in genes.keys():
		genes[line[0]]=[]
		genes[line[0]].append((1-(1-float(line[2]))**.5)*(numpermute+1)-1)
	else:
		genes[line[0]].append((1-(1-float(line[2]))**.5)*(numpermute+1)-1)
outfile = open("Heiko_700_SeedScores_merge",'w')
for gene in genes:
	outfile.write(gene+"\t"+str((sum(genes[gene])+1)/(numpermute*len(genes[gene])+1))+"\n")
outfile.close()


outfile = open("Heiko_400_NetStats_merge",'w')

diseaseresults = open("heiko_input_400_P62_NetStats",'r').readlines()
x = float(diseaseresults[1].split("\t")[1])
os.system("cat heiko_input_400_P*_permutedDirectEdgeCount > heiko_input_400_permutedDirectEdgeCount")
random=[]
file = open("heiko_input_400_permutedDirectEdgeCount",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Direct Edge Count\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[2].split("\t")[2])
os.system("cat heiko_input_400_P*_permutedSeedDirectDegreesMean > heiko_input_400_permutedSeedDirectDegreesMean")
random=[]
file = open("heiko_input_400_permutedSeedDirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Direct Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[3].split("\t")[2])
os.system("cat heiko_input_400_P*_permutedIndirectDegreesMean > heiko_input_400_permutedIndirectDegreesMean")
random=[]
file = open("heiko_input_400_permutedIndirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Indirect Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[4].split("\t")[2])
os.system("cat heiko_input_400_P*_permutedCIdegreesMeanPermut > heiko_input_400_permutedCIdegreesMeanPermut")
random=[]
file = open("heiko_input_400_permutedCIdegreesMeanPermut",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("CI Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

outfile.close()

os.system("cat heiko_input_400_P*_SeedScores > heiko_input_400_SeedScores")
file = open("heiko_input_400_SeedScores",'r')
genes = {}
for line in file:
	line = line.strip("\n").split()
	if line[0]=="GENE" or line[2]=="NA":
		continue
	if not line[0] in genes.keys():
		genes[line[0]]=[]
		genes[line[0]].append((1-(1-float(line[2]))**.5)*(numpermute+1)-1)
	else:
		genes[line[0]].append((1-(1-float(line[2]))**.5)*(numpermute+1)-1)
outfile = open("Heiko_400_SeedScores_merge",'w')
for gene in genes:
	outfile.write(gene+"\t"+str((sum(genes[gene])+1)/(numpermute*len(genes[gene])+1))+"\n")
outfile.close()

sys.exit()

outfile = open("Heiko_full_NetStats_merge",'w')

diseaseresults = open("heiko_input_full_P62_NetStats",'r').readlines()
x = float(diseaseresults[1].split("\t")[1])
os.system("cat heiko_input_full_P*_permutedDirectEdgeCount > heiko_input_full_permutedDirectEdgeCount")
random=[]
file = open("heiko_input_full_permutedDirectEdgeCount",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Direct Edge Count\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[2].split("\t")[2])
os.system("cat heiko_input_full_P*_permutedSeedDirectDegreesMean > heiko_input_full_permutedSeedDirectDegreesMean")
random=[]
file = open("heiko_input_full_permutedSeedDirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Direct Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[3].split("\t")[2])
os.system("cat heiko_input_full_P*_permutedIndirectDegreesMean > heiko_input_full_permutedIndirectDegreesMean")
random=[]
file = open("heiko_input_full_permutedIndirectDegreesMean",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("Seed Indirect Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

x = float(diseaseresults[4].split("\t")[2])
os.system("cat heiko_input_full_P*_permutedCIdegreesMeanPermut > heiko_input_full_permutedCIdegreesMeanPermut")
random=[]
file = open("heiko_input_full_permutedCIdegreesMeanPermut",'r')
for line in file:
	random.append(float(line.strip("\n")))
outfile.write("CI Degrees Mean\t"+str(float(len([i for i in random if i>x])+1)/float(len(random)+1))+"\n")

outfile.close()
