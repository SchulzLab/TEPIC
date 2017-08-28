#This script combines TF affinities and Peak features with gene expression data (that is present in a tab delimited format)
#Note that both files are expected to have a header
#arg1: TEPIC gene-TF scores
#arg2: Gene expression estimates. assumed format: GeneID \tab Value
#arg3: Output file
#--geneIDs: Column of gene IDs in the file specified with arg2
#--expressionC Column of gene expression estimates in the file specified with arg2
#--filterIDs File containing geneIDs that should be used for intersection (File with one column only, no header)

import sys
import os
import argparse

def main():
	parser=argparse.ArgumentParser(prog="annotateTSS.py")
	parser.add_argument("affinities",nargs=1,help="TEPIC gene-TF scores")
	parser.add_argument("expression",nargs=1,help="Gene expression file")
	parser.add_argument("output",nargs=1,help="File to write the combined data to")
	parser.add_argument("--geneIDs",nargs="?",help="Position of the gene IDs in the expression file",default=0,type=int)
	parser.add_argument("--expressionC",nargs="?",help="Position of the gene expression estimate in the expression file",default=1,type=int)
	parser.add_argument("--filterIDs",nargs="?",help="File containing gene IDs that should be considered",default=None)
	args=parser.parse_args()
	
	print("Reading TF affinities from file: "+args.affinities[0])
	tfFile=open(args.affinities[0],"r")
	tfFileheader=tfFile.readline().strip()
	tfDict={}
	tfKeys=set()
	for l in tfFile:
		s=l.split()
		tfDict[s[0]]=l.strip()
		tfKeys.add(s[0])
	tfFile.close()

	print("Reading Gene expression from file: "+args.expression[0])
	if (args.geneIDs != 0):
		print("Gene ID are retrieved from column "+str(args.geneIDs))
	if (args.expressionC != 1):
		print("Gene expression estimates are retrieved from column "+str(args.expressionC))
	expFile=open(args.expression[0],"r")
	expFileheader=expFile.readline().strip()
	expDict={}
	expKeys=set()
	for l in expFile:
		s=l.split()
		if ("." in s[0]):
			expDict[s[0].split(".")[0]]=str(s[args.expressionC])
			expKeys.add(s[0].split(".")[0])	
		else:
			expDict[s[0]]=str(s[args.expressionC])
			expKeys.add(s[0])
	expFile.close()
	keys=expKeys.intersection(tfKeys)

	if (args.filterIDs !=None):
		filterSet=set()
		print("Loading IDs that should be used for filtering")
		filterFile=open(args.filterIDs,"r")
		for l in filterFile:
			if ("." in l):
				filterSet.add(l.split(".")[0])
			else:
				filterSet.add(l.strip())			
		filterFile.close()
		keys=keys.intersection(filterSet)
	
	print("Overlapping gene IDs: "+str(len(keys)))
	print("Writing integrated data to "+str(args.output[0]))
	outfile=open(args.output[0],"w")
	outfile.write(tfFileheader+"\tExpression\n")
	for key in keys:
		outfile.write(tfDict[key]+"\t"+expDict[key]+"\n")
	outfile.close()

main()	
