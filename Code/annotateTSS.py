import sys
import string
import operator
import math
import argparse
from decimal import Decimal

#Computing per gene TF affinities

#Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
def readGTF(filename):
	gtf=open(sys.argv[1],"r")
	tss={}
	for l in gtf:
		s=l.split()
		if (len(s) >=9):
			if (s[2]=="gene"):	
				if (s[6]=="+"):
					tss[s[9]]=(s[0].replace("chr",""),int(s[3]))
				else:
					tss[s[9]]=(s[0].replace("chr",""),int(s[4]))
	gtf.close()
	return tss

#Reads the txt file containing TF-scores. Extracts the regions of open chromatin.
#They are returned as a dictionary(key: #chromosom, item:[(start,end)])
def readOC_Region(filename):
	tfpa=open(filename,"r")
	oC={}
	for l in tfpa:
		s=l.split()[0]
		ds=s.split(":")
		if (len(ds)>=2):
			se=ds[1].split("-")
			if (not oC.has_key(ds[0].replace("chr",""))):
				oC[ds[0].replace("chr","")]=[(int(se[0]),int(se[1]))]			
			else:
				oC[ds[0].replace("chr","")]+=[(int(se[0]),int(se[1]))]
	tfpa.close()
	return oC


def aggregateAffinity(old,new,factor):
	for i in xrange(0,len(old)-1):
		old[i]=float(old[i])+factor*float(new[i])
	return old


def extractTF_Affinity(openRegions,genesInOpenChromatin,filename,genePositions,openChromatin,expDecay):
	geneAffinities={}
	tfpa=open(filename,"r")
	tfpa.readline()
	for l in tfpa:
		s=l.split()
		middles=s[0].split(":")[1].split("-")
		middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
		if (genesInOpenChromatin.has_key(s[0])):
			for geneID in genesInOpenChromatin[s[0]]:
				tss=genePositions[geneID][1]
				factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
				if(s[0] in openRegions):
					if (expDecay):
						if (geneAffinities.has_key(geneID)):
							geneAffinities[geneID]=aggregateAffinity(geneAffinities[geneID],s[1:],factor)
						else:
							numbers=s[1:]
							for i in xrange(0,len(numbers)-1):
								numbers[i]=float(factor)*float(numbers[i])
							geneAffinities[geneID]=numbers				
					else:
						if (geneAffinities.has_key(geneID)):
							geneAffinities[geneID]=aggregateAffinity(geneAffinities[geneID],s[1:],1.0)
						else:
							geneAffinities[geneID]=s[1:]

	tfpa.close()
	return geneAffinities


def tfIndex(filename):
	tfpa=open(filename,"r")
	l=tfpa.readline()
	tfpa.close()
	return l.split()

def createAffinityFile(affinities,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	output.write(header+'\n')
	for Gene in tss.keys():
		line=""
		if (affinities.has_key(Gene)):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in xrange(0,len(tfNames)):
				line+='\t'+str(0.0)
		output.write(line+'\n')
	output.close()


def createSparseFile(affinities,tfNames,filename,tss,number):
	if (len(tfNames) < number):
		number=len(tfNames)
		print("Warning: The value of sparseRep is to large, representation will contain all possible TFs")
	output=open(filename,"w")
	header="geneID\tTF\tAffinity\n"
	output.write(header)
	for Gene in tss.keys():
		tfList=[]
		if (affinities.has_key(Gene)):
			geneID=str(Gene.replace("\"","").replace(";","").split(".")[0])
			temp=affinities[Gene]
			for i in xrange(0,len(tfNames)):
				tfList=tfList+[((tfNames[i],float(temp[i])))]
			tfList.sort(key=lambda tup:tup[1],reverse=True)
			for i in xrange(0,number):	
				output.write(str(geneID)+"\t"+str(tfList[i][0])+"\t"+str(tfList[i][1])+"\n")
	output.close()


def makeTupels(values,names):
	l=[]
	for i in xrange(0,len(values)-1):
		l+=[(names[i],values[i])]
	return l

def main():
	parser=argparse.ArgumentParser(prog="annotateTSSV2.py")
	parser.add_argument("gtf",nargs=1,help="Genome annotation file")
	parser.add_argument("affinity",nargs=1,help="TRAP generated TF Affinity file")
	parser.add_argument("--geneViewAffinity",nargs="?",help="Name of the gene view affinity files. If this is not specified, the prefix of the input files will be used.",default="")
	parser.add_argument("--windows",nargs="?",help="Size of the considered window around the TSS. Default is 3000.",default=3000,type=int)
	parser.add_argument("--decay",nargs="?",help="True if exponential decay should be used, False otherwise. Default is True",default="True")
	parser.add_argument("--signalScale",nargs="?",help="If the name of the scaled affinity file is provied, a Gene view file is computed for those Affinity values.",default="")
	parser.add_argument("--sparseRep",nargs="?",help="Number of top TFs that should be contained in the sparse representation",default=0,type=int)
	args=parser.parse_args() 

	prefixs=args.affinity[0].split(".")
	prefix=prefixs[0]
	if (args.geneViewAffinity==""):
		args.geneViewAffinity=prefix+"_Affinity_Gene_View.txt"

	if (args.decay.upper()=="FALSE") or (args.decay=="0"):
		decay=False
	else:
		decay=True

	#Check arguments
	
	#Extract TSS of GTF files
	tss=readGTF(args.gtf[0])
	#Load open chromatin positions from TF-Affinity file
	oC=readOC_Region(args.affinity[0])
	#Create a TF name index
	tfNames=tfIndex(args.affinity[0])
	shift=int(args.windows/2)
	#Determine gene windows in open chromatin regions
	genesInOpenChromatin={}
	usedRegions=set()
	for gene in tss.keys():
		if (oC.has_key(tss[gene][0])):
			for tupel in oC[tss[gene][0]]:
				#Right border of window <= Right border of open chromatin
				if (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift >= tupel[0]):
					#Left border of window >= Left border of open chromatin ==> Window inside open chromatin
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))					
				#Right border of window >= Left border of open chromatin ==> Window enters open chromatin in the 3' end and stops in the tss window
				elif (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift < tupel[0]) and (tss[gene][1]+shift > tupel[0]):
					
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))
				#Right border of window > Right border of open chromatin
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift < tupel[0]):
					#Left border of window <= Left border of open chromatin ==> Window is larger than open chromatin
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))
				#Left border of window <= Right border of open chromain ==> Window enters open chromatin in the 5' end stops in the tss window
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift >= tupel[0]) and (tss[gene][1]-shift < tupel[1]):
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))
	

	#Extract bound transcription factors
	affinities=extractTF_Affinity(usedRegions,genesInOpenChromatin,args.affinity[0],tss,oC,decay)
	if (decay):
		createAffinityFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Affinity_Gene_View.txt"),tss)	
		if (args.sparseRep != 0):
			createSparseFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Sparse_Affinity_Gene_View.txt"),tss,args.sparseRep)
	else:
		createAffinityFile(affinities,tfNames,args.geneViewAffinity,tss)
		if (args.sparseRep != 0):
			createSparseFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Sparse_Affinity_Gene_View.txt"),tss,args.sparseRep)



	scaledAffinities={}
	if (args.signalScale != ""):
		scaledAffinities=extractTF_Affinity(usedRegions,genesInOpenChromatin,args.signalScale,tss,oC,decay)
		if (decay):
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Affinity_Gene_View.txt"),tss)
			if (args.sparseRep != 0):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Sparse_Affinity_Gene_View.txt"),tss,args.sparseRep)
		else:
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Affinity_Gene_View.txt"),tss)	
			if (args.sparseRep != 0):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Sparse_Scaled_Affinity_Gene_View.txt"),tss,args.sparseRep)





main()
