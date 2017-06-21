#This script combines TF affinities and Peak features with gene expression data (that is present in a tab delimited two column format)
#Note that both files are expected to have a header

import sys
import os

def main():
	
	if (len(sys.argv) < 4):
		print("Usage: python integrateData.py <TF-Affinities> <Gene Expression> <Integrated Data>")
		exit()
	print("Reading TF affinities from file: "+sys.argv[1])
	tfFile=open(sys.argv[1],"r")
	tfFileheader=tfFile.readline().strip()
	tfDict={}
	tfKeys=set()
	for l in tfFile:
		s=l.split()
		tfDict[s[0]]=l.strip()
		tfKeys.add(s[0])
	tfFile.close()

	print("Reading Gene expression from file: "+sys.argv[2])
	expFile=open(sys.argv[2],"r")
	expFileheader=expFile.readline().strip()
	expDict={}
	expKeys=set()
	for l in expFile:
		s=l.split()
		expDict[s[0]]=str(s[1])
		expKeys.add(s[0])
	expFile.close()
	keys=expKeys.intersection(tfKeys)
	
	print("Overlapping gene IDs: "+str(len(keys)))
	print("Writing integrated data")
	outfile=open(sys.argv[3],"w")
	outfile.write(tfFileheader+"\tExpression\n")
	for key in keys:
		outfile.write(tfDict[key]+"\t"+expDict[key]+"\n")
	outfile.close()

main()	
