import sys
import argparse
import os

def readFile(filename,geneSet,affinityDict):
	temp=set()
	infile=open(filename,"r")
	header=infile.readline()
	for l1 in infile:
		s1=l1.split()
		if (affinityDict.has_key(s1[0])):
			for i in xrange(1,len(s1)):
				affinityDict[s1[0]][i-1]=float(affinityDict[s1[0]][i-1])+float(s1[i])
		else:		
			affinityDict[s1[0]]=s1[1:]
		temp.add(s1[0])
		geneSet.add(s1[0])
	print("Found "+str(len(temp))+" genes")
	infile.close()
	return (header,geneSet,affinityDict)

def divideGeneExpression(affinityDict,gene,counter):
	temp=gene
	s=affinityDict[gene]
	for i in xrange(0,len(s)):
		temp+="\t"+str(float(s[i])/counter)
	return temp+"\n"

def checkHeader(headers):
	for i in xrange(0,len(headers)-2):
		if (headers[i]!=headers[i+1]):
			return False
	return True
		
def main():
	parser=argparse.ArgumentParser(prog="computeMeanRatioTFAffinities.py")
	parser.add_argument("group1",nargs=1,help="Path to the folder containing TF affinities calculated for group 1")
	parser.add_argument("group2",nargs=1,help="Path to the folder containing TF affinities calculated for group 2")
	parser.add_argument("ogroup1",nargs=1,help="Name of the file to store the mean affinities of group 1")
	parser.add_argument("ogroup2",nargs=1,help="Name of the file to store the mean affinities of group 2")
	parser.add_argument("oratios",nargs=1,help="Name of the file to store the affinity ratios between groups 1 and 2")
	parser.add_argument("scaled",nargs=1,help="Flag indiciating wheter scaled affinity features should be used", default="False")
	parser.add_argument("peakFeatures",nargs=1,help="Flag indiciating wheter peakFeatures should be used", default="False")
	args=parser.parse_args()
	genes=set()
	affinityDictG1={}
	counterG1=0
	affinityDictG2={}
	counterG2=0
	headers=[]
	if((args.peakFeatures[0].upper()=="TRUE") | (args.peakFeatures[0].upper()=="T")):
		if (args.scaled[0]=="True"):
			suffix="Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt"
		else:
			suffix="Peak_Features_Affinity_Gene_View_Filtered.txt"
	else:
		if (args.scaled[0]=="True"):
			suffix="Signal_Feature_Affinity_Gene_View_Filtered.txt"
		else:
			suffix="Affinity_Gene_View_Filtered.txt"
	for group1file in os.listdir(args.group1[0]):
		if (suffix in group1file):
			print("Processing file: "+group1file)
			(header,genes,affinityDictG1)=readFile(args.group1[0]+"/"+group1file,genes,affinityDictG1)
			headers+=[header]
			counterG1+=1
	for group2file in os.listdir(args.group2[0]):
		if (suffix in group2file):
			print("Processing file: "+group2file)
			(header,genes,affinityDictG2)=readFile(args.group2[0]+"/"+group2file,genes,affinityDictG2)
			headers+=[header]
			counterG2+=1
	for gene in genes:
		if (affinityDictG1.has_key(gene)):
			continue
		else:
			affinityDictG1[gene]=[0.0]*(len(headers[0].split())-1)
	for gene in genes:
		if (affinityDictG2.has_key(gene)):
			continue
		else:
			affinityDictG2[gene]=[0.0]*(len(headers[0].split())-1)

	if (len(genes)==0):
		print("There are no genes available") 
		exit()
	print("Total number of genes: "+str(len(genes)))

	if (checkHeader(headers)):
		print("Calculating mean affinities for group 1")
		outfile1=open(args.ogroup1[0],"w")
		outfile1.write(str(headers[0]))
		for geneID in genes:
			outfile1.write(divideGeneExpression(affinityDictG1,geneID,counterG1))
		outfile1.close()
		print("Calculating mean affinities for group 2")
		outfile2=open(args.ogroup2[0],"w")
		outfile2.write(str(headers[0]))
		for geneID in genes:
			outfile2.write(divideGeneExpression(affinityDictG2,geneID,counterG2))
		outfile2.close()
	else:
		print("Error,processed files do not contain data for the same TFs")

	print("Processing Means of group 1")
	infile1=open(args.ogroup1[0],"r")
	header=infile1.readline().strip()
	s1dict={}
	for l1 in infile1:
		s1=l1.split()
		s1dict[s1[0]]=s1[1:]
	infile1.close()

	print("Processing Means of group 2")
	infile2=open(args.ogroup2[0],"r")
	header=infile2.readline().strip()
	s2dict={}
	for l2 in infile2:
		s2=l2.split()
		s2dict[s2[0]]=s2[1:]
	infile2.close()
	pc=1.0

	print("Computing affinity ratios")
	outfile3=open(args.oratios[0],"w")
	outfile3.write(str(headers[0]))
	for key in genes:
		s1=s1dict[key]
		s2=s2dict[key]
		temp=str(key)
		for i in xrange(0,len(s1)):
			temp+="\t"+str(((float(s1[i])+pc)/(float(s2[i])+pc)))
		outfile3.write(temp+"\n")
	outfile3.close()

main()
