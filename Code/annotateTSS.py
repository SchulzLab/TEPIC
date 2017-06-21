import sys
import string
import operator
import math
import argparse
import random
from decimal import Decimal

#Computing per gene TF affinities
#Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
def readGTF(filename):
	gtf=open(sys.argv[1],"r")
	open(filename,"r")
	identifier="start_codon"
	for l in gtf:
		s=l.split()
		if (len(s) >=9):
			if (s[2]=="gene"):
				identifier="gene"
				break
	gtf.close()
	gtf=open(filename,"r")
	tss={}
	for l in gtf:
		s=l.split()
		if (len(s) >=9):
			if (s[2]==identifier):
				if (s[6]=="+"):
					if (tss.has_key(s[9])):
						if (int(s[3]) < tss[s[9]][1]):
							tss[s[9]]=(s[0].replace("chr",""),(int(s[3]),int(s[4])))
					else:
						tss[s[9]]=(s[0].replace("chr",""),(int(s[3]),int(s[4])))
				else:
					if (tss.has_key(s[9])):
						if (int(s[4]) > tss[s[9]][1]):
							tss[s[9]]=(s[0].replace("chr",""),(int(s[4]),int(s[3])))
					else:
						tss[s[9]]=(s[0].replace("chr",""),(int(s[4]),int(s[3])))
	gtf.close()
	return tss,identifier

#Reads the txt file containing TF-scores. Extracts the regions of open chromatin.
#They are returned as a dictionary(key: #chromosom, item:[(start,end)])
def readOC_Region(filename):
	tfpa=open(filename,"r")
	tfpa.readline()
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

def extractTF_Affinity(openRegions,genesInOpenChromatin,filename,genePositions,openChromatin,expDecay,geneBody,peakFeatures,lengthNormalisation,motifLength):
	geneAffinities={}
	numberOfPeaks={}
	totalPeakLength={}
	tfpa=open(filename,"r")
	tfpa.readline()
	if (not geneBody):
		for l in tfpa:
			s=l.split()
			middles=s[0].split(":")[1].split("-")
			middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
			length=float(int(middles[1])-int(middles[0]))
			if (genesInOpenChromatin.has_key(s[0])):
				for geneID in genesInOpenChromatin[s[0]]:
					if(s[0] in openRegions):
						tss=genePositions[geneID][1][0]
						if (expDecay):
							factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
							if (geneAffinities.has_key(geneID)):
								if (lengthNormalisation):	
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
								else:
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
								totalPeakLength[geneID]+=length*factor
								numberOfPeaks[geneID]+=factor
							else:
								numbers=map(lambda x: float(x)*float(factor),s[1:])
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.div,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
								else:
									geneAffinities[geneID]=numbers
								totalPeakLength[geneID]=length*factor
								numberOfPeaks[geneID]=factor
						else:
							if (geneAffinities.has_key(geneID)):
								numberOfPeaks[geneID]+=1.0
								totalPeakLength[geneID]+=length
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
								else:
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
							else:
								numberOfPeaks[geneID]=1.0
								totalPeakLength[geneID]=length
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
								else:
									geneAffinities[geneID]=map(lambda x: float(x),s[1:])
	else:
		for l in tfpa:
			s=l.split()
			middles=s[0].split(":")[1].split("-")
			middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
			length=float(int(middles[1])-int(middles[0]))
			if (genesInOpenChromatin.has_key(s[0])):
				for geneID in genesInOpenChromatin[s[0]]:
					if(s[0] in openRegions):
						tss=genePositions[geneID][1][0]
						tts=genePositions[geneID][1][1]
						if (tss < tts):
							if (middle < tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneAffinities.has_key(geneID)):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
										numberOfPeaks[geneID]+=factor
										totalPeakLength[geneID]+=(factor*length)
									else:
										numbers=map(lambda x:float(factor)*float(x),s[1:])
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.div,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=numbers				
										numberOfPeaks[geneID]=factor
										totalPeakLength[geneID]=length*factor
								else:
									if (geneAffinities.has_key(geneID)):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
										numberOfPeaks[geneID]+=1.0
										totalPeakLength[geneID]+=length
									else:
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=map(lambda x: float(x),s[1:])
										numberOfPeaks[geneID]=1.0
										totalPeakLength[geneID]=length
							else:
								if (geneAffinities.has_key(geneID)):
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
									else:
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x), s[1:]))
									totalPeakLength[geneID]+=length
									numberOfPeaks[geneID]+=1.0
								else:
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
									else:
										geneAffinities[geneID]=map(lambda x: flota(x),s[1:])
									totalPeakLength[geneID]=length
									numberOfPeaks[geneID]=1.0
						else:
							if (middle > tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneAffinities.has_key(geneID)):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
										numberOfPeaks[geneID]+=factor
										totalPeakLength[geneID]+=(factor*length)
									else:
										numbers=map(lambda x: float(x)*float(factor),s[1:])
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.div,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=numbers				
										numberOfPeaks[geneID]=factor
										totalPeakLength[geneID]=length*factor
								else:
									if (geneAffinities.has_key(geneID)):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
										numberOfPeaks[geneID]+=1.0
										totalPeakLength[geneID]+=length
									else:
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=map(lambda x: float(x), s[1:])
										numberOfPeaks[geneID]=1.0
										totalPeakLength[geneID]=length
							else:
								if (geneAffinities.has_key(geneID)):
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
									else:
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
									numberOfPeaks[geneID]+=1.0
									totalPeakLength[geneID]+=length
								else:
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.div,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
									else:
											geneAffinities[geneID]=map(lambda x: float(x), s[1:])
									numberOfPeaks[geneID]=1.0
									totalPeakLength[geneID]=length
	tfpa.close()
	return geneAffinities,numberOfPeaks,totalPeakLength

def generate_Peak_Coverage_Features(openRegions,genesInOpenChromatin,filename,genePositions,openChromatin,expDecay,geneBody):
	perBaseDNaseSignal={}
	tfpa=open(filename,"r")
	tfpa.readline()
	if (not geneBody):
		for l in tfpa:
			s=l.split()
			peakPos=s[0]+":"+s[1]+"-"+s[2]
			middle=int(((float(s[2])-float(s[1]))/2)+float(s[1]))
			if (genesInOpenChromatin.has_key(peakPos)):
				for geneID in genesInOpenChromatin[peakPos]:
					if(peakPos in openRegions):
						tss=genePositions[geneID][1][0]
						if (expDecay):
							factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
							if (perBaseDNaseSignal.has_key(geneID)):
								perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
							else:
								perBaseDNaseSignal[geneID]=float(s[3])*factor
						else:
							if (perBaseDNaseSignal.has_key(geneID)):
								perBaseDNaseSignal[geneID]+=float(s[3])
							else:
								perBaseDNaseSignal[geneID]=float(s[3])	
	else:
		for l in tfpa:
			s=l.split()
			peakPos=s[0]+":"+s[1]+"-"+s[2]
			middle=int(((float(s[2])-float(s[1]))/2)+float(s[1]))
			length=float(int(s[2])-int(s[1]))
			if (genesInOpenChromatin.has_key(s[0])):
				for geneID in genesInOpenChromatin[s[0]]:
					if(s[0] in openRegions):
						tss=genePositions[geneID][1][0]
						tts=genePositions[geneID][1][1]
						if (tss < tts):
							if (middle < tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (perBaseDNaseSignal.has_key(geneID)):
										perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
									else:
										perBaseDNaseSignal[geneID]=float(s[3])*factor
								else:
									if (perBaseDNaseSignal.has_key(geneID)):
										perBaseDNaseSignal[geneID]+=float(s[3])
									else:
										perBaseDNaseSignal[geneID]=float(s[3])
							else:
								if (perBaseDNaseSignal.has_key(geneID)):
									perBaseDNaseSignal[geneID]+=float(s[3])
								else:
									perBaseDNaseSignal[geneID]=float(s[3])
						else:
							if (middle > tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (perBaseDNaseSignal.has_key(geneID)):
										perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
									else:
										perBaseDNaseSignal[geneID]=float(s[3])*factor
								else:
									if (perBaseDNaseSignal.has_key(geneID)):
										perBaseDNaseSignal[geneID]+=float(s[3])
									else:
										perBaseDNaseSignal[geneID]=float(s[3])
							else:
								if (perBaseDNaseSignal.has_key(geneID)):
									perBaseDNaseSignal[geneID]+=float(s[3])
								else:
									perBaseDNaseSignal[geneID]=float(s[3])

	tfpa.close()
	return perBaseDNaseSignal


def tfIndex(filename):
	tfpa=open(filename,"r")
	l=tfpa.readline()
	tfpa.close()
	return l.split()

#Creates an affinity file that contains only TF affinities, no additional features
def createAffinityFileAffintiesOnly(affinities,tfNames,filename,tss):
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

#Creates an affinity file that contains affinities, peak counts and peak length information
def createAffinityFileAffinitiesPeakCountsLength(affinities,peakCounts,peakLength,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	header+="\tPeak_Counts\tPeak_Length"
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
		if (peakCounts.has_key(Gene)):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (peakLength.has_key(Gene)):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

#Creates an affinity file that contains affinities, peak counts, peak length, and peak signal information
def createAffinityFileAffinitiesPeakCountsLengthSignal(affinities,peakCounts,peakLength,peakSignal,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	header+="\tPeak_Counts\tPeak_Length\tPeak_Signal"
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
		if (peakCounts.has_key(Gene)):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (peakLength.has_key(Gene)):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		if (peakSignal.has_key(Gene)):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

#Creates an affinity file that contains affinities, and peak signal information
def createAffinityFileAffinitiesSignal(affinities,peakSignal,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	header+="\tPeak_Signal"
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
		if (peakSignal.has_key(Gene)):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()


def createSparseFile(affinities,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID\tTF\tAffinity\n"
	output.write(header)
	for Gene in tss.keys():
		if (affinities.has_key(Gene)):
			geneID=str(Gene.replace("\"","").replace(";","").split(".")[0])
			temp=affinities[Gene]
			for i in xrange(0,len(tfNames)):
				if (float(temp[i]) > 0):
					output.write(str(geneID)+"\t"+str(tfNames[i])+"\t1\n")
	output.close()


def makeTupels(values,names):
	l=[]
	for i in xrange(0,len(values)-1):
		l+=[(names[i],values[i])]
	return l

def generate_Motif_Length(affinityFile,motifFile):
	motifDict={}
	motifList=[]
	affinityF=open(affinityFile,"r")
	header=affinityF.readline().split()
	affinityF.close()
	if (motifFile!=None):
		motifF=open(motifFile,"r")
		for l in motifF:
			s=l.split()
			motifDict[s[0]]=int(s[1])
		motifF.close()
	for tf in header:
		if (motifDict.has_key(tf.upper())):
			motifList+=[motifDict[tf.upper()]]
		else:
			motifList+=[0]
	return motifList

def main():
	parser=argparse.ArgumentParser(prog="annotateTSS.py")
	parser.add_argument("gtf",nargs=1,help="Genome annotation file")
	parser.add_argument("affinity",nargs=1,help="TRAP generated TF Affinity file")
	parser.add_argument("--geneViewAffinity",nargs="?",help="Name of the gene view affinity files. If this is not specified, the prefix of the input files will be used.",default=None)
	parser.add_argument("--windows",nargs="?",help="Size of the considered window around the TSS. Default is 3000.",default=3000,type=int)
	parser.add_argument("--decay",nargs="?",help="True if exponential decay should be used, False otherwise. Default is True",default="True")
	parser.add_argument("--signalScale",nargs="?",help="If the name of the scaled affinity file is provied, a Gene view file is computed for those Affinity values.",default=None)
	parser.add_argument("--sparseRep",nargs="?",help="Flag to be set if a sparse representation should be generated. This should only be used with a filtered set of TF affinities",default="False")
	parser.add_argument("--geneBody",nargs="?",help="True if the entire gene body should be screened for TF binding",default="False")
	parser.add_argument("--peakCoverage",nargs="?",help="File containing the per base DNase1 signal in the peaks",default=None)
	parser.add_argument("--additionalPeakFeatures",nargs="?",help="True if additional features based on peak count and peak length should be computed, Default is False",default="False")	
	parser.add_argument("--randomizePerGene",nargs="?",help="Flag indicating whether, in addition to the standard output, a matrix with randomized features per gene should be generated. Default is False.",default="False")
	parser.add_argument("--normaliseLength",nargs="?",help="Normalises the TF affinities with the total peak length. Default is False.",default="False")
	parser.add_argument("--motifLength",nargs="?",help="File containing the length of the used motifs. Used to adapt the length normalisation such that long motifs are not downweighted compared to short ones",default=None)
	args=parser.parse_args() 

	prefixs=args.affinity[0].split(".")
	prefix=prefixs[0]
	if (args.geneViewAffinity==None):
		args.geneViewAffinity=prefix+"_Affinity_Gene_View.txt"

	if (args.decay.upper()=="FALSE") or (args.decay=="0") or (args.decay.upper()=="F"):
		decay=False
	else:
		decay=True

	if (args.geneBody.upper()=="FALSE") or (args.geneBody=="0") or (args.geneBody.upper()=="F"):
		geneBody=False
	else:
		geneBody=True

	if (args.additionalPeakFeatures.upper()=="FALSE") or (args.additionalPeakFeatures=="0") or (args.additionalPeakFeatures.upper()=="F"):
		addPeakF=False
		addPeakFT=False
	else:
		addPeakF=True
		addPeakFT=True

	if (args.randomizePerGene.upper()=="FALSE") or (args.randomizePerGene=="0") or (args.randomizePerGene.upper()=="F"):
		randomizePerGene=False
	else:
		randomizePerGene=True

	if (args.normaliseLength.upper()=="FALSE") or (args.normaliseLength.upper()=="0") or (args.normaliseLength.upper()=="F"):
		normaliseLength=False
	else:
		normaliseLength=True
		addPeakFT=True
	
	if (args.sparseRep.upper()=="FALSE") or (args.sparseRep=="0") or (args.sparseRep.upper()=="F"):
		sparseRep=False
	else:
		sparseRep=True
	
	#Check arguments
	#Extract TSS of GTF files
	tss,identifier=readGTF(args.gtf[0])                                                                                                                                   
	if (identifier=="start_codon"):
		geneBody=False
		print("Gene Body parameter forced to be false, due to incompatible gene annotation")
	#Generate List of motif lengths
	motifLengths=generate_Motif_Length(args.affinity[0],args.motifLength)	
	#Load open chromatin positions from TF-Affinity file
	oC=readOC_Region(args.affinity[0])
	#Create a TF name index
	tfNames=tfIndex(args.affinity[0])
	shift=int(args.windows/2)
	#Determine gene windows in open chromatin regions
	genesInOpenChromatin={}
	usedRegions=set()
	for gene in tss.keys():
		#Define window borders here		
		if (not geneBody):
			leftBorder=tss[gene][1][0]-shift
			rightBorder=tss[gene][1][0]+shift
		else:
			if (int(tss[gene][1][0]) < int(tss[gene][1][1])):
				leftBorder=tss[gene][1][0]-shift
				rightBorder=tss[gene][1][1]
			else:
				leftBorder=tss[gene][1][1]
				rightBorder=tss[gene][1][0]+shift
		if (oC.has_key(tss[gene][0])):
			for tupel in oC[tss[gene][0]]:
				#Right border of window <= Right border of open chromatin
				if (rightBorder <= tupel[1]) and (leftBorder >= tupel[0]):
					#Left border of window >= Left border of open chromatin ==> Window inside open chromatin
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))					
				#Right border of window >= Left border of open chromatin ==> Window enters open chromatin in the 5' end and stops in the tss window
				elif (rightBorder <= tupel[1]) and (leftBorder < tupel[0]) and (rightBorder > tupel[0]):
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))
				#Right border of window > Right border of open chromatin
				elif (rightBorder > tupel[1]) and (leftBorder < tupel[0]):
					#Left border of window <= Left border of open chromatin ==> Window is larger than open chromatin
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))
				#Left border of window <= Right border of open chromain ==> Window enters open chromatin in the 3' end stops in the tss window
				elif (rightBorder > tupel[1]) and (leftBorder >= tupel[0]) and (leftBorder < tupel[1]):
					if (genesInOpenChromatin.has_key(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))):		
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]+=[gene]
					else:
						genesInOpenChromatin[tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])]=[gene]
					usedRegions.add(tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1]))

	#Extract bound transcription factors
	affinities,numberOfPeaks,peakLength=extractTF_Affinity(usedRegions,genesInOpenChromatin,args.affinity[0],tss,oC,decay,geneBody,addPeakFT,normaliseLength,motifLengths)
	#Generate Peak based features
	if (args.peakCoverage != None):
		perBaseCoverage=generate_Peak_Coverage_Features(usedRegions,genesInOpenChromatin,args.peakCoverage,tss,oC,decay,geneBody)

	#Generate Output
	if (decay):
		if (addPeakF):
			createAffinityFileAffinitiesPeakCountsLength(affinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Peak_Features_Affinity_Gene_View.txt"),tss)	
		else:
			createAffinityFileAffintiesOnly(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Affinity_Gene_View.txt"),tss)		
		if (sparseRep):
			createSparseFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Sparse_Affinity_Gene_View.txt"),tss)
	else:
		if (addPeakF):
			createAffinityFileAffinitiesPeakCountsLength(affinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Peak_Features_Affinity_Gene_View.txt"),tss)
		else:
			createAffinityFileAffintiesOnly(affinities,tfNames,args.geneViewAffinity,tss)
		if (sparseRep):
			createSparseFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Sparse_Affinity_Gene_View.txt"),tss)

	scaledAffinities={}
	if (args.signalScale != None):
		scaledAffinities,numberOfPeaks,peakLength=extractTF_Affinity(usedRegions,genesInOpenChromatin,args.signalScale,tss,oC,decay,geneBody,addPeakF,normaliseLength,motifLengths)
		if (decay):
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLength(scaledAffinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Peak_Features_Affinity_Gene_View.txt"),tss)
			else:
				createAffinityFileAffintiesOnly(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Affinity_Gene_View.txt"),tss)
			if (sparseRep):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Sparse_Affinity_Gene_View.txt"),tss)
		else:
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLength(scaledAffinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Peak_Features_Affinity_Gene_View.txt"),tss)	
			else:
				createAffinityFileAffintiesOnly(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Affinity_Gene_View.txt"),tss)	
			if (sparseRep):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Sparse_Scaled_Affinity_Gene_View.txt"),tss)
	
	if (args.peakCoverage != None):
		if (decay):
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLengthSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt"),tss)	
			else:
				createAffinityFileAffinitiesSignal(affinities,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Signal_Feature_Affinity_Gene_View.txt"),tss)	
		else:
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLengthSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Three_Peak_Based_Features_Affinity_Gene_View.txt"),tss)
			else:
				createAffinityFileAffinitiesSignal(affinities,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Signal_Feature_Affinity_Gene_View.txt"),tss)

main()
