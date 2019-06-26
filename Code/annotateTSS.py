import sys
import string
import operator
from operator import itemgetter, add
import math
import argparse
import random
from decimal import Decimal
from SortedCollection import SortedCollection
from pprint import pprint

#Computing per gene TF affinities
#Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
def readGTF(filename,tA):
	gtf=open(sys.argv[1],"r")
	open(filename,"r")
	identifier="start_codon"
	idIndex=9
	if (tA==True):
		identifier="transcript"
		idIndex=11
	else:
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
		if (len(s)>=idIndex):
			if (s[2]==identifier):
				if (s[6]=="+"):
					if (s[idIndex] in tss):
						if (int(s[3]) < tss[s[idIndex]][1]):
							tss[s[idIndex]]=(s[0].replace("chr",""),(int(s[3]),int(s[4])))
					else:
						tss[s[idIndex]]=(s[0].replace("chr",""),(int(s[3]),int(s[4])))
				else:
					if (s[idIndex] in tss):
						if (int(s[4]) > tss[s[idIndex]][1]):
							tss[s[idIndex]]=(s[0].replace("chr",""),(int(s[4]),int(s[3])))
					else:
						tss[s[idIndex]]=(s[0].replace("chr",""),(int(s[4]),int(s[3])))
	gtf.close()
	return tss,identifier

#Reads the txt file containing TF-scores. Extracts the regions of open chromatin.
#They are returned as a dictionary(key: #chromosom, item:[(start,end)])
def readOC_Region(filename):
	tfpa=open(filename,"r")
	tfpa.readline()
	oC={}
	counter=1
	for l in tfpa:
		s=l.split()[0]
		ds=s.split(":")
		if (len(ds)>=2):
			chrom=ds[0].replace("chr","")
			se=ds[1].split("-")
			if chrom not in oC:
				oC[chrom]=SortedCollection(key=itemgetter(1))
			oC[chrom].insert_right((counter,int(se[0]),int(se[1])))
			counter+=1
	tfpa.close()
	return oC


		
#Extract the TF affinities and compute peak length and count features from the TF affinity file
def extractTF_Affinity(openRegions,genesInOpenChromatin,filename,genePositions,openChromatin,expDecay,geneBody,peakFeatures,lengthNormalisation,motifLength):
	geneAffinities={}
	numberOfPeaks={}
	totalPeakLength={}
	tfpa=open(filename,"r")
	headerlength=len(tfpa.readline().split())
	if (not geneBody):
		for l in tfpa:
			s=l.split()
			if (headerlength != len(s) -1):
				print("Header of TF affinity file is not correctly formatted. Number of columns does not match the header")
				exit()
			middles=s[0].split(":")[1].split("-")
			middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
			length=float(int(middles[1])-int(middles[0]))
			if (s[0] in genesInOpenChromatin):
				for geneID in genesInOpenChromatin[s[0]]:
					if(s[0] in openRegions):
						tss=genePositions[geneID][1][0]
						if (expDecay):
							factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
							if (geneID in geneAffinities):
								if (lengthNormalisation):	
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
								else:
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
								totalPeakLength[geneID]+=length*factor
								numberOfPeaks[geneID]+=factor
							else:
								numbers=map(lambda x: float(x)*float(factor),s[1:])
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.truediv,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
								else:
									geneAffinities[geneID]=numbers
								totalPeakLength[geneID]=length*factor
								numberOfPeaks[geneID]=factor
						else:
							if (geneID in geneAffinities):
								numberOfPeaks[geneID]+=1.0
								totalPeakLength[geneID]+=length
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
								else:
									geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
							else:
								numberOfPeaks[geneID]=1.0
								totalPeakLength[geneID]=length
								if (lengthNormalisation):
									geneAffinities[geneID]=map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
								else:
									geneAffinities[geneID]=map(lambda x: float(x),s[1:])
	else:
		for l in tfpa:
			s=l.split()
			if (headerlength != len(s) -1):
				print("Header of TF affinity file is not correctly formatted. Number of columns does not match the header")
				exit()

			middles=s[0].split(":")[1].split("-")
			middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
			length=float(int(middles[1])-int(middles[0]))
			if (s[0] in genesInOpenChromatin):
				for geneID in genesInOpenChromatin[s[0]]:
					if(s[0] in openRegions):
						tss=genePositions[geneID][1][0]
						tts=genePositions[geneID][1][1]
						if (tss < tts):
							if (middle < tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneID in geneAffinities):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
										numberOfPeaks[geneID]+=factor
										totalPeakLength[geneID]+=(factor*length)
									else:
										numbers=map(lambda x:float(factor)*float(x),s[1:])
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.truediv,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=numbers				
										numberOfPeaks[geneID]=factor
										totalPeakLength[geneID]=length*factor
								else:
									if (geneID in geneAffinities):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
										numberOfPeaks[geneID]+=1.0
										totalPeakLength[geneID]+=length
									else:
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=map(lambda x: float(x),s[1:])
										numberOfPeaks[geneID]=1.0
										totalPeakLength[geneID]=length
							else:
								if (geneID in geneAffinities):
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
									else:
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x), s[1:]))
									totalPeakLength[geneID]+=length
									numberOfPeaks[geneID]+=1.0
								else:
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
									else:
										geneAffinities[geneID]=map(lambda x: float(x),s[1:])
									totalPeakLength[geneID]=length
									numberOfPeaks[geneID]=1.0
						else:
							if (middle > tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneID in geneAffinities):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(lambda x: factor*float(x),s[1:]),map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: factor*float(x),s[1:]))
										numberOfPeaks[geneID]+=factor
										totalPeakLength[geneID]+=(factor*length)
									else:
										numbers=map(lambda x: float(x)*float(factor),s[1:])
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.truediv,numbers,map(lambda x: factor*(length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=numbers				
										numberOfPeaks[geneID]=factor
										totalPeakLength[geneID]=length*factor
								else:
									if (geneID in geneAffinities):
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
										else:
											geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
										numberOfPeaks[geneID]+=1.0
										totalPeakLength[geneID]+=length
									else:
										if (lengthNormalisation):
											geneAffinities[geneID]=map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
										else:
											geneAffinities[geneID]=map(lambda x: float(x), s[1:])
										numberOfPeaks[geneID]=1.0
										totalPeakLength[geneID]=length
							else:
								if (geneID in geneAffinities):
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))) 
									else:
										geneAffinities[geneID]=map(operator.add,geneAffinities[geneID],map(lambda x: float(x),s[1:]))
									numberOfPeaks[geneID]+=1.0
									totalPeakLength[geneID]+=length
								else:
									if (lengthNormalisation):
										geneAffinities[geneID]=map(operator.truediv,map(float,s[1:]),map(lambda x: (length-x+1) if (length-x+1 > 0) else 1, motifLength))
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
			if (peakPos in genesInOpenChromatin):
				for geneID in genesInOpenChromatin[peakPos]:
					if(peakPos in openRegions):
						tss=genePositions[geneID][1][0]
						if (expDecay):
							factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
							if (geneID in perBaseDNaseSignal):
								perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
							else:
								perBaseDNaseSignal[geneID]=float(s[3])*factor
						else:
							if (geneID in perBaseDNaseSignal):
								perBaseDNaseSignal[geneID]+=float(s[3])
							else:
								perBaseDNaseSignal[geneID]=float(s[3])	
	else:
		for l in tfpa:
			s=l.split()
			peakPos=s[0]+":"+s[1]+"-"+s[2]
			middle=int(((float(s[2])-float(s[1]))/2)+float(s[1]))
			length=float(int(s[2])-int(s[1]))
			if (peakPos in genesInOpenChromatin):
				for geneID in genesInOpenChromatin[peakPos]:
					if(peakPos in openRegions):
						tss=genePositions[geneID][1][0]
						tts=genePositions[geneID][1][1]
						if (tss < tts):
							if (middle < tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneID in perBaseDNaseSignal):
										perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
									else:
										perBaseDNaseSignal[geneID]=float(s[3])*factor
								else:
									if (geneID in perBaseDNaseSignal):
										perBaseDNaseSignal[geneID]+=float(s[3])
									else:
										perBaseDNaseSignal[geneID]=float(s[3])
							else:
								if (geneID in perBaseDNaseSignal):
									perBaseDNaseSignal[geneID]+=float(s[3])
								else:
									perBaseDNaseSignal[geneID]=float(s[3])
						else:
							if (middle > tss):
								if (expDecay):
									factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
									if (geneID in perBaseDNaseSignal):
										perBaseDNaseSignal[geneID]+=(float(s[3])*factor)
									else:
										perBaseDNaseSignal[geneID]=float(s[3])*factor
								else:
									if (geneID in perBaseDNaseSignal):
										perBaseDNaseSignal[geneID]+=float(s[3])
									else:
										perBaseDNaseSignal[geneID]=float(s[3])
							else:
								if (geneID in perBaseDNaseSignal):
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
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
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
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)
		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
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
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)
		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignal):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()


def createConformationDataOutputPeaks(peakCounts,peakLength,peakSignal,peakCountsLR,peakLengthLR,peakSignalLR,fNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	header+="\tPeak_Counts\tPeak_Length\tPeak_Signal\tPeak_CountsLR\tPeak_LengthLR\tPeak_SignalLR"
	output.write(header+'\n')
	for Gene in tss.keys():
		line=str(Gene.replace("\"","").replace(";","").split(".")[0])
		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignal):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		if (Gene in peakCountsLR):
			line+='\t'+str(peakCountsLR[Gene])
		else:
			line+='\t0'
		if (Gene in peakLengthLR):
			line+='\t'+str(peakLengthLR[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignalLR):
			line+='\t'+str(peakSignalLR[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

def createConformationDataOutputPeaksCombined(peakCounts,peakLength,peakSignal,peakCountsLR,peakLengthLR,peakSignalLR,fNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	header+="\tPeak_Counts\tPeak_Length\tPeak_Signal"
	output.write(header+'\n')
	for Gene in tss.keys():
		line=str(Gene.replace("\"","").replace(";","").split(".")[0])
		counts=0.0
		if (Gene in peakCounts):
			counts=float(peakCounts[Gene])
		if (Gene in peakCountsLR):
			counts+=float(peakCountsLR[Gene])
		line+='\t'+str(counts)

		length=0.0
		if (Gene in peakLength):
			length=float(peakLength[Gene])
		if (Gene in peakLengthLR):
			length+=float(peakLengthLR[Gene])
		line+='\t'+str(length)

		signal=0.0
		if (Gene in peakSignal):
			signal=float(peakSignal[Gene])
		if (Gene in peakSignalLR):
			signal+=float(peakSignalLR[Gene])
		line+='\t'+str(signal)

		output.write(line+'\n')
	output.close()


def createConformationDataOutputAffinitiesPeaks(affinities,affinitiesLR,peakCounts,peakLength,peakSignal,peakCountsLR,peakLengthLR,peakSignalLR,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	for element in tfNames:
		header+='\t'+"LR_"+str(element)
	header+="\tPeak_Counts\tPeak_Length\tPeak_Signal\tPeak_CountsLR\tPeak_LengthLR\tPeak_SignalLR"
	output.write(header+'\n')
	for Gene in tss.keys():
		line=""
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)

		if (Gene in affinitiesLR):
			for entry in affinitiesLR[Gene]:
				line+='\t'+str(entry)
		else:
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)

		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignal):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'

		if (Gene in peakCountsLR):
			line+='\t'+str(peakCountsLR[Gene])
		else:
			line+='\t0'
		if (Gene in peakLengthLR):
			line+='\t'+str(peakLengthLR[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignalLR):
			line+='\t'+str(peakSignalLR[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()


def createConformationDataOutputAffinities(affinities,affinitiesLR,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID"
	for element in tfNames:
		header+='\t'+str(element)
	for element in tfNames:
		header+='\t'+"LR_"+str(element)
	output.write(header+'\n')
	for Gene in tss.keys():
		line=""
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)

		if (Gene in affinitiesLR):
			for entry in affinitiesLR[Gene]:
				line+='\t'+str(entry)
		else:
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)
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
		if (Gene in affinities):
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for entry in affinities[Gene]:
				line+='\t'+str(entry)
		else:
			line=str(Gene.replace("\"","").replace(";","").split(".")[0])
			for i in range(0,len(tfNames)):
				line+='\t'+str(0.0)
		if (Gene in peakSignal):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

def createPeakScoreFileExcludingSignal(affinities,peakCounts,peakLength,filename,tss):
	output=open(filename,"w")
	header="geneID"
	header+="\tPeak_Counts\tPeak_Length"
	output.write(header+'\n')
	for Gene in tss.keys():
		line=str(Gene.replace("\"","").replace(";","").split(".")[0])
		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

def createPeakScoreFileIncludingSignal(affinities,peakCounts,peakLength,peakSignal,filename,tss):
	output=open(filename,"w")
	header="geneID"
	header+="\tPeak_Counts\tPeak_Length\tPeak_Signal"
	output.write(header+'\n')
	for Gene in tss.keys():
		line=str(Gene.replace("\"","").replace(";","").split(".")[0])
		if (Gene in peakCounts):
			line+='\t'+str(peakCounts[Gene])
		else:
			line+='\t0'
		if (Gene in peakLength):
			line+='\t'+str(peakLength[Gene])
		else:
			line+='\t0'
		if (Gene in peakSignal):
			line+='\t'+str(peakSignal[Gene])
		else:
			line+='\t0'
		output.write(line+'\n')
	output.close()

#Creates a sparse representation file containing binary gene TF-relations
def createSparseFile(affinities,tfNames,filename,tss):
	output=open(filename,"w")
	header="geneID\tTF\tAffinity\n"
	output.write(header)
	for Gene in tss.keys():
		if (Gene in  affinities):
			geneID=str(Gene.replace("\"","").replace(";","").split(".")[0])
			temp=affinities[Gene]
			for i in range(0,len(tfNames)):
				if (float(temp[i]) > 0):
					output.write(str(geneID)+"\t"+str(tfNames[i])+"\t1\n")
	output.close()


def makeTupels(values,names):
	l=[]
	for i in range(0,len(values)-1):
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
		if (tf.upper() in motifDict):
			motifList+=[motifDict[tf.upper()]]
		else:
			motifList+=[0]
	return motifList


#Reads the DNA contact file
#Returns two sorted collections as dictionaries(key: #chromsom, item(l1s,l1e,l2s,l2e) and #chromsom, item(l2s,l2e,l1s,l1e).
def read_Intra_Loop_Regions_Collection(filename):
	left_region_collection = {}
	right_region_collection = {}
	with open(filename) as hi_c_file:
		loopID=0
		for line in hi_c_file:
			line = line.split()
			if (line[0] == line[3]):
				loopID += 1
				chromosome = line[0].replace("chr","")
				if chromosome not in left_region_collection:
					left_region_collection[chromosome] = SortedCollection(key=itemgetter(1))
				if chromosome not in right_region_collection:
					right_region_collection[chromosome] = SortedCollection(key=itemgetter(1))	
				left_region_collection[chromosome].insert_right((loopID, int(line[1]), int(line[2]), int(line[4]), int(line[5])))
				right_region_collection[chromosome].insert_right((loopID, int(line[4]), int(line[5]), int(line[1]), int(line[2])))
	return left_region_collection, right_region_collection


#Intersect peaks with Loops
def getGenesInLongRangeWindows(annotations,regions_collection,shift):
	gene_regions = {}
	for geneID in annotations.keys():
		gene_regions[geneID]=SortedCollection(key=itemgetter(1))
		left_Border=annotations[geneID][1][0]-shift
		right_Border=annotations[geneID][1][0]+shift
		chromosome=annotations[geneID][0]
		if (chromosome in regions_collection):
			selectedLoops = regions_collection[chromosome]
			try:                                                                                                     
				left_item = selectedLoops.find_lt(left_Border)
			except ValueError:
				try:
					left_item = selectedLoops.find_ge(left_Border)
				except ValueError:
					left_item = None
			else:
				if left_item[2] < left_Border:
					try:
						left_item = selectedLoops.find_ge(left_Border)
					except ValueError:
						left_item = None
			try:
				right_item = selectedLoops.find_le(right_Border)
			except ValueError:
				right_item = None
	           # Check if target interval is valid
			if left_item is not None and right_item is not None:
				left_index = selectedLoops.index(left_item)
				right_index = selectedLoops.index(right_item)
				if left_index <= right_index:
	          # Copy regions in target interval
					for i in xrange(left_index, right_index + 1):
						gene_regions[geneID].insert_right(selectedLoops[i])
	return gene_regions	
		

def retrieveOverlappingOpenSites(dhs_collection, loops, annotation):
	filtered_dhs = {}
	for geneID in loops.keys():
		chromosome=annotation[geneID][0]
		if (chromosome in dhs_collection):	
			regions = dhs_collection[chromosome]
			filtered_dhs[geneID] = get_intersecting_regions(regions,loops[geneID])
	return filtered_dhs


def get_intersecting_regions(a_regions, b_collection):
	intersection_a = SortedCollection(key=itemgetter(1))
	for a_region in a_regions:
		try:
			left_boundary = b_collection.find_lt_index(a_region[1])
		except ValueError:
			try:
				left_boundary = b_collection.find_ge_index(a_region[1])
			except ValueError:
				left_boundary = len(b_collection)
		else:
			if b_collection[left_boundary][2] < a_region[1]:
				left_boundary += 1

		curr_index = left_boundary
		if curr_index < len(b_collection) and b_collection[curr_index][1] <= a_region[2]:
			intersection_a.insert_right(a_region)
			curr_index += 1

	b_collection.key=itemgetter(3)
	for a_region in a_regions:
		try:
			left_boundary = b_collection.find_lt_index(a_region[1])
		except ValueError:
			try:
				left_boundary = b_collection.find_ge_index(a_region[1])
			except ValueError:
				left_boundary = len(b_collection)
		else:
			if b_collection[left_boundary][4] < a_region[1]:
				left_boundary += 1

		curr_index = left_boundary
		if curr_index < len(b_collection) and b_collection[curr_index][3] <= a_region[2]:
			intersection_a.insert_right(a_region)
			curr_index += 1
	b_collection.key=itemgetter(1)
	return intersection_a


def merge_gene_regions(gene_regions, add_gene_regions,annotation):
	for geneID in annotation.keys():
		if geneID not in gene_regions.keys():
			gene_regions[geneID] = SortedCollection(key=itemgetter(1))
		if geneID in add_gene_regions.keys():
			for region in add_gene_regions[geneID]:
				try:
					gene_regions[geneID].find(region[1])
				except ValueError:
					gene_regions[geneID].insert_right(region)
	return gene_regions

def remove_gene_regions(gene_regions, promoterDHSlist,annotation):
	for geneID in gene_regions.keys():
		loopDHS=gene_regions[geneID]
		for region in list(loopDHS):
			identifier=str(annotation[geneID][0])+":"+str(region[1])+"-"+str(region[2])
			if (identifier in promoterDHSlist):
				if (geneID in promoterDHSlist[identifier]):
					try:
						gene_regions[geneID].remove(region)
					except ValueError:
						pass
	return gene_regions


#convert it to a set of all used DHS (chr:start-end) and to a dict of geneIDs-> [DHS as chr:start-end]
def convert_Conformation_Data(filteredConformation_DHSs, annotation):
	usedLongRangeDHS=set()
	usedLongRangeGeneIDs={}
	for geneID in filteredConformation_DHSs.keys():
		if (geneID in filteredConformation_DHSs):
			chrom=annotation[geneID][0]
			temp=list(filteredConformation_DHSs[geneID])
			for element in temp:
				identifier=str(chrom)+":"+str(element[1])+"-"+str(element[2])
				if identifier not in usedLongRangeGeneIDs:
					usedLongRangeGeneIDs[identifier]=[]
				usedLongRangeGeneIDs[identifier]+=[geneID]
				usedLongRangeDHS.add(identifier)
	return usedLongRangeDHS,usedLongRangeGeneIDs

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
	parser.add_argument("--onlyPeakFeatures",nargs="?",help="Generates an additional output file that contains only peak based features per gene. Default is False.",default="False")
	parser.add_argument("--transcript",nargs="?",help="Extract the position of transcripts from the gtf file and generate a transcript based annotation instead of a gene centric one. Default is False.",default="False")
	parser.add_argument("--lwindows",nargs="?",help="Size of the considered loop window around the TSS. Default is 5000.",default=5000,type=int)
	parser.add_argument("--conformationData",nargs="?",help="File holding chromatin contact information",default=None)
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
	
	if (args.sparseRep.upper()=="FALSE") or (args.sparseRep.upper()=="0") or (args.sparseRep.upper()=="F"):
		sparseRep=False
	else:
		sparseRep=True

	if (args.onlyPeakFeatures.upper()=="FALSE") or (args.onlyPeakFeatures.upper()=="0") or (args.onlyPeakFeatures.upper()=="F"):
		onlyPeakFeatures=False
	else:
		onlyPeakFeatures=True

	if (args.transcript.upper()=="FALSE") or (args.transcript.upper()=="0") or (args.transcript.upper()=="F"):
		transcriptAnnotation=False
	else:
		transcriptAnnotation=True
	
	#Check arguments
	#Extract TSS of GTF files
	tss,identifier=readGTF(args.gtf[0],transcriptAnnotation)   
	if (identifier=="start_codon"):
		geneBody=False
		print("Gene Body parameter forced to be false, due to incompatible gene annotation")
	#Generate List of motif lengths
	motifLengths=generate_Motif_Length(args.affinity[0],args.motifLength)	
	#Load open chromatin positions from TF-Affinity file
	oC=readOC_Region(args.affinity[0])
	if (not oC):
		print("No TF affinities provided in "+args.affinity[0]+".")
		return 0
	#Create a TF name index
	tfNames=tfIndex(args.affinity[0])
	shift=int(args.windows/2)

	###Generate promoter only features###
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
		chrom=tss[gene][0]
		if (chrom in oC):
			try:
				left_item = oC[chrom].find_lt(leftBorder)
			except ValueError:
				try:
					left_item = oC[chrom].find_ge(leftBorder)
				except ValueError:
					left_item = None
			else:
				if left_item[2] < leftBorder:
					try:
						left_item = oC[chrom].find_ge(leftBorder)
					except ValueError:
						left_item = None
			try:
				right_item = oC[chrom].find_le(rightBorder)
			except ValueError:
				right_item = None
			if left_item is not None and right_item is not None:
				left_index = oC[chrom].index(left_item)
				right_index = oC[chrom].index(right_item)
				if left_index <= right_index:
					for i in range(left_index, right_index + 1):
						identifier=str(chrom)+":"+str(oC[chrom][i][1])+"-"+str(oC[chrom][i][2])
						if identifier in genesInOpenChromatin:
							genesInOpenChromatin[identifier]+=[gene]
						else:
							genesInOpenChromatin[identifier]=[gene]
						usedRegions.add(str(chrom)+":"+str(oC[chrom][i][1])+"-"+str(oC[chrom][i][2]))

	#Extract bound transcription factors
	affinities,numberOfPeaks,peakLength=extractTF_Affinity(usedRegions,genesInOpenChromatin,args.affinity[0],tss,oC,decay,geneBody,addPeakFT,normaliseLength,motifLengths)

	#Generate Peak based features
	if (args.peakCoverage != None):
		perBaseCoverage=generate_Peak_Coverage_Features(usedRegions,genesInOpenChromatin,args.peakCoverage,tss,oC,decay,geneBody)


	###Generate 3D based features###
	if (args.conformationData != None):
		(regions_left_collection, regions_right_collection) = read_Intra_Loop_Regions_Collection(args.conformationData)

		leftGenes=getGenesInLongRangeWindows(tss,regions_left_collection,float(args.lwindows)/2.0)
		rightGenes=getGenesInLongRangeWindows(tss,regions_right_collection,float(args.lwindows)/2.0)

		leftOverlap=retrieveOverlappingOpenSites(oC, leftGenes, tss)
		rightOverlap=retrieveOverlappingOpenSites(oC, rightGenes, tss)
		totalConformation_list=merge_gene_regions(leftOverlap, rightOverlap,tss)
		filteredConformation_DHSs=remove_gene_regions(totalConformation_list, genesInOpenChromatin,tss)
		#convert it to a set of all used DHS (chr:start-end) and to a dict of geneIDs-> [DHS as chr:start-end]
		usedRegionsLongRange,genesInLongRange=convert_Conformation_Data(filteredConformation_DHSs,tss)
	
		affinitiesLR,numberOfPeaksLR,peakLengthLR=extractTF_Affinity(usedRegionsLongRange,genesInLongRange,args.affinity[0],tss,oC,False,False,addPeakFT,normaliseLength,motifLengths)
		if (args.peakCoverage != None):
			perBaseCoverageLR=generate_Peak_Coverage_Features(usedRegionsLongRange,genesInLongRange,args.peakCoverage,tss,oC,False,geneBody)
	

	#Generate Output
	if (decay):
		if (addPeakF):
			createAffinityFileAffinitiesPeakCountsLength(affinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Peak_Features_Affinity_Gene_View.txt"),tss)
			if (onlyPeakFeatures):
				createPeakScoreFileExcludingSignal(affinities,numberOfPeaks,peakLength,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Peak_Features_Only_Gene_View.txt"),tss)
		else:
			createAffinityFileAffintiesOnly(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Affinity_Gene_View.txt"),tss)		
		if (sparseRep):
			createSparseFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Sparse_Affinity_Gene_View.txt"),tss)
	else:
		if (addPeakF):
			createAffinityFileAffinitiesPeakCountsLength(affinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Peak_Features_Affinity_Gene_View.txt"),tss)
			if (onlyPeakFeatures):
				createPeakScoreFileExcludingSignal(affinities,numberOfPeaks,peakLength,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Peak_Features_Only_Gene_View.txt"),tss)
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
				if (onlyPeakFeatures):
					createPeakScoreFileExcludingSignal(scaledAffinities,numberOfPeaks,peakLength,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Peak_Features_Only_Gene_View.txt"),tss)
			else:
				createAffinityFileAffintiesOnly(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Affinity_Gene_View.txt"),tss)
			if (sparseRep):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Sparse_Affinity_Gene_View.txt"),tss)
		else:
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLength(scaledAffinities,numberOfPeaks,peakLength,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Peak_Features_Affinity_Gene_View.txt"),tss)
				if (onlyPeakFeatures):
					createPeakScoreFileExcludingSignal(scaledaffinities,numberOfPeaks,peakLength,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Peak_Features_Only_Gene_View.txt"),tss)
			else:
				createAffinityFileAffintiesOnly(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Affinity_Gene_View.txt"),tss)	
			if (sparseRep):
				createSparseFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Sparse_Scaled_Affinity_Gene_View.txt"),tss)
	
	if (args.peakCoverage != None):
		if (decay):
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLengthSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt"),tss)
				if (onlyPeakFeatures):
					createPeakScoreFileIncludingSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Three_Peak_Based_Features_Only_Gene_View.txt"),tss)
			else:
				createAffinityFileAffinitiesSignal(affinities,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Signal_Feature_Affinity_Gene_View.txt"),tss)	
		else:
			if (addPeakF):
				createAffinityFileAffinitiesPeakCountsLengthSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Three_Peak_Based_Features_Affinity_Gene_View.txt"),tss)
				if (onlyPeakFeatures):
					createPeakScoreFileIncludingSignal(affinities,numberOfPeaks,peakLength,perBaseCoverage,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Three_Peak_Based_Features_Only_Gene_View.txt"),tss)
			else:
				createAffinityFileAffinitiesSignal(affinities,perBaseCoverage,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Signal_Feature_Affinity_Gene_View.txt"),tss)

	if (args.conformationData != None):
		if (addPeakF):
			if (args.peakCoverage != None):
				if (onlyPeakFeatures):
					createConformationDataOutputPeaks(numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLR,peakLengthLR,perBaseCoverageLR,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Three_Peak_Based_Features_Double_Gene_View.txt"),tss)
				else:
					createConformationDataOutputPeaks(numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLR,peakLengthLR,perBaseCoverageLR,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Three_Peak_Based_Features_Double_Gene_View.txt"),tss)
					createConformationDataOutputAffinitiesPeaks(affinities,affinitiesLR,numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLR,peakLengthLR,perBaseCoverageLR,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Affinity_Three_Peak_Based_Features_Gene_View.txt"),tss)
			elif (args.signalScale != None):
				affinitiesLRS,numberOfPeaksLRS,peakLengthLRS=extractTF_Affinity(usedRegionsLongRange,genesInLongRange,args.signalScale,tss,oC,False,False,addPeakFT,normaliseLength,motifLengths)
				if (onlyPeakFeatures):
					createConformationDataOutputPeaks(numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLRS,peakLengthLRS,perBaseCoverageLRS,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Three_Peak_Based_Features_Double_Gene_View.txt"),tss)
				else:
					createConformationDataOutputPeaks(numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLRS,peakLengthLRS,perBaseCoverageLRS,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Three_Peak_Based_Features_Double_Gene_View.txt"),tss)
					createConformationDataOutputAffinitiesPeaks(affinities,affinitiesLRS,numberOfPeaks,peakLength,perBaseCoverage,numberOfPeaksLRS,peakLengthLRS,perBaseCoverageLRS,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Affinity_Three_Peak_Based_Features_Gene_View.txt"),tss)
		else:
			createConformationDataOutputAffinities(affinities,affinitiesLR,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Conformation_Data_Affinity_Gene_View.txt"),tss)

main()
