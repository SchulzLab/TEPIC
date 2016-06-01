import sys
import string
import operator
import math
import argparse
import utils
from decimal import Decimal

#Computing per gene TF affinities

#Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
def readGTF(filename):
	gtf=open(filename,"r")
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


def readTFPA(filename):
	f = open(filename,"r")
	f.readline()
	tfpa = {}
	for l in f:
		s=l.split()
		tfpa[s[0]] = s
	f.close()
	return tfpa


def aggregateAffinity(old,new,factor):
	for i in xrange(0,len(old)-1):
		old[i]=float(old[i])+factor*float(new[i])
	return old


def extractTF_Affinity(genesInOpenChromatin,filename,genePositions,expDecay, loopsactivated, looptable):
	geneAffinities={}
	tfpa=readTFPA(filename)
	for posKey in tfpa:
		s = tfpa[posKey]
		middles=s[0].split(":")[1].split("-")
		middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
		if (genesInOpenChromatin.has_key(s[0])):
			for geneID in genesInOpenChromatin[s[0]][0]:
				tss=genePositions[geneID][1]
				if (expDecay):
					factor=math.exp(-(float(float(abs(tss-middle))/5000.0)))
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
				if(loopsactivated):	#maybe change this here
					loopid = genesInOpenChromatin[s[0]][1]
					loop = looptable[loopid]
					loopsite = genesInOpenChromatin[s[0]][2]
					regions = []
					if(loopsite):# we are on the left loop-site, now also add the regions from the counterpart loop-site
						regions = loop[1]
					else:
						regions = loop[0]
					for region in countersiteocregions:
						newpos = s[0].split(":")[0]+region[0]+"-"+region[1]
						snew = tfpa[newpos]
						geneAffinities[geneID] = aggregateAffinity(geneAffinities[geneID],snew[1:],1.0)

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

def makeTupels(values,names):
	l=[]
	for i in xrange(0,len(values)-1):
		l+=[(names[i],values[i])]
	return l
	
	
def intersectRegions(oc, loops):
	intersectedOC = {}
	looptable = {}
	for chrKey in oc:
		if(not intersectedOC.hasKey(chrKey)):
			intersectedOC[chrKey] = []
		for octupel in oc[chrKey]:
			if(loops.hasKey(chrKey)):
				chrLoops = loops[chrKey]
				for loop in chrLoops:
					# open chromatin region in...
					#	[  ----  ]
					if(octupel[0] >= loop[1] and octupel[1] <= loop[2])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], True) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# --[------  ]
					elif(octupel[0] <= loop[1] and octupel[1] <= loop[2] && octupel[1] >= loop[1])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], True) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# 	[  ------]--	
					elif(octupel[0] >= loop[1] and octupel[0] <= loop[2] && octupel[1] >= loop[2])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], True) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# --[--------]--	
					elif(octupel[0] < loop[1] and octupel[1] > loop[2])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], True) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# same procedure for second, counterpart loop-site
					#	[  ----  ]
					elif(octupel[0] >= loop[3] and octupel[1] <= loop[4])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], False) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# --[------  ]
					elif(octupel[0] <= loop[3] and octupel[1] <= loop[4] && octupel[1] >= loop[1])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], False) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# 	[  ------]--	
					elif(octupel[0] >= loop[3] and octupel[0] <= loop[4] && octupel[1] >= loop[2])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], False) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# --[--------]--	
					elif(octupel[0] < loop[3] and octupel[1] > loop[4])):
						intersectedOC[chrKey] += (octupel[0], octupel[1], loop[0], False) # also append the loop-id
						if(not looptable.hasKey(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
						
	return (intersectedOC, looptable)


def filterLoops(loops, resolution):
	for chrKey in loops:
		chrLoops = loops[chrKey]
		for loop in chrLoops:
			if((loop[2]-loop[1]) != resolution):
				chrLoops.remove(loop)
	
				
def findLoopsNearbyGenes(tss, loops, loopwindows, usemiddle):
	geneLoops = {}
	for geneID in tss:
		startpos = tss[geneID][1]
		chrLoops = loops[tss[geneID][0]]
		if(not geneLoops.has_key(geneID)):
			geneLoops[geneID] = []
		for loop in chrLoops:
			foundleftone = False
			foundrightone = False
			if(usemiddle):
				
				leftmiddle = int(((float(loop[2])-float(loop[1]))/2)+float(loop[1]))
				if(abs(middle - startpos) <= loopwindows):
					foundleftone = True
					
				rightmiddle = int(((float(loop[3])-float(loop[4]))/2)+float(loop[3]))
				if (abs(middle - startpos) <= loopwindows):
					foundrightone = True
				
				if(foundleftone and foundrightone):
					if(leftmiddle <= rightmiddle):
						geneLoops[geneID].append((loop[0], True))
					
				elif(foundleftone):
					geneLoops[geneID].append((loop[0], True))
				
				elif(foundrightone):
					geneLoops[geneID].append((loop[0], False))
					
			else:
				leftpos = 0
				if((abs(loop[1] - startpos) <= loopwindows) or (abs(loop[2] - startpos) <= loopwindows)):
					foundleftone = True
					if(abs(loop[1] - startpos) <= abs(loop[2] - startpos)):
						leftpos = abs(loop[1] - startpos)
					else:
						leftpos = abs(loop[2] - startpos)
				
				rightpos = 0
				if((abs(loop[3] - startpos) <= loopwindows) or (abs(loop[4] - startpos) <= loopwindows)):
					foundrightone = True
					if(abs(loop[3] - startpos) <= abs(loop[4] - startpos)):
						rightpos = abs(loop[3] - startpos)
					else:
						rightpos = abs(loop[4] - startpos)
					
				if(foundleftone and foundrightone):
					if (leftpos <= rightpos):
						geneLoops[geneID].append((loop[0], True))
				elif(foundleftone):
					geneLoops[geneID].append((loop[0], True))
				elif(foundrightone):
					geneLoops[geneID].append((loop[0], False))
				
	return geneLoops
					

def main():
	parser=argparse.ArgumentParser(prog="annotateTSSV2.py")
	parser.add_argument("gtf",nargs=1,help="Genome annotation file")
	parser.add_argument("affinity",nargs=1,help="TRAP generated TF Affinity file")
	parser.add_argument("--geneViewAffinity",nargs="?",help="Name of the gene view affinity files. If this is not specified, the prefix of the input files will be used.",default="")
	parser.add_argument("--windows",nargs="?",help="Size of the considered window around the TSS. Default is 3000.",default=3000,type=int)
	parser.add_argument("--decay",nargs="?",help="True if exponential decay should be used, False otherwise. Default is True",default="True")
	parser.add_argument("--signalScale",nargs="?",help="If the name of the scaled affinity file is provied, a Gene view file is computed for those Affinity values.",default="")
	parser.add_argument("--loopfile",nargs="?",help="If the name of the loop file is provied, all open chromatin regions will be intersected with loop regions around the TSS of each gene.",default="")
	parser.add_argument("--loopwindows",nargs="?",help="Defines the window-size around the TSS in which all loops are considered for intersecting with openChromatin regions.",default=50000,type=int)
	parser.add_argument("--resolution",nargs="?",help="Defines the Hi-C resolution of the loops which should be considered. Uses the smallest one found if the a resolution is not available in the loopfile.",default=5000,type=int)
	parser.add_argument("--usemiddle",nargs="?",help="Defines whether to use the middle of a loop to decide if a loop lies withing a window or the edges.",default="False")
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
	loopsactivated = False
	looptable = {}
	#Extract loops of Hi-C loopfile
	if(args.loopfile != ""):
		loopwindows = args.loopwindows
		resolution = args.resolution
		loopsactivated = True
		
		#looptable = utils.readIntraLoopsLookupTable(args.loopfile) # maybe use this later on
		loops = utils.readIntraLoops(args.loopfile)
		#filter loops and keep user-defined resolution only
		filterLoops(loops, resolution)
		#intersect openchromatin regions with all loopregions in TSS windows
		intersectResults = intersectRegions(oC, loops)
		oC = intersectresults[0]
		looptable = intersectresults[1]
		
	
	#Determine gene windows in open chromatin regions
	genesInOpenChromatin={}
	usedRegions=set()
	for gene in tss.keys():
		if (oC.has_key(tss[gene][0])):
			for tupel in oC[tss[gene][0]]:
				loci = tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])
				#Right border of window <= Right border of open chromatin
				if (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift >= tupel[0]):
					#Left border of window >= Left border of open chromatin ==> Window inside open chromatin
					if (genesInOpenChromatin.has_key(loci)):		
						genesInOpenChromatin[loci][0]+=[gene]
					else:
						if(loopsactivated):
							genesInOpenChromatin[loci]=([gene],tupel[2],tupel[3])
						else:
							genesInOpenChromatin[loci]=([gene])
					usedRegions.add(loci)					
				#Right border of window >= Left border of open chromatin ==> Window enters open chromatin in the 3' end and stops in the tss window
				elif (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift < tupel[0]) and (tss[gene][1]+shift > tupel[0]):
					
					if (genesInOpenChromatin.has_key(loci)):		
						genesInOpenChromatin[loci][0]+=[gene]
					else:
						if(loopsactivated):
							genesInOpenChromatin[loci]=([gene],tupel[2],tupel[3])
						else:
							genesInOpenChromatin[loci]=([gene])
						
					usedRegions.add(loci)
				#Right border of window > Right border of open chromatin
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift < tupel[0]):
					#Left border of window <= Left border of open chromatin ==> Window is larger than open chromatin
					if (genesInOpenChromatin.has_key(loci)):		
						genesInOpenChromatin[loci][0]+=[gene]
					else:
						if(loopsactivated):
							genesInOpenChromatin[loci]=([gene],tupel[2],tupel[3])
						else:
							genesInOpenChromatin[loci]=([gene])
					usedRegions.add(loci)
				#Left border of window <= Right border of open chromain ==> Window enters open chromatin in the 5' end stops in the tss window
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift >= tupel[0]) and (tss[gene][1]-shift < tupel[1]):
					if (genesInOpenChromatin.has_key(loci)):		
						genesInOpenChromatin[loci][0]+=[gene]
					else:
						if(loopsactivated):
							genesInOpenChromatin[loci]=([gene],tupel[2],tupel[3])
						else:
							genesInOpenChromatin[loci]=([gene])
					usedRegions.add(loci)
	

	#Extract bound transcription factors
	affinities=extractTF_Affinity(genesInOpenChromatin,args.affinity[0],tss,decay,loopsactivated, looptable)
	if (decay):
		createAffinityFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Affinity_Gene_View.txt"),tss)	
	else:
		createAffinityFile(affinities,tfNames,args.geneViewAffinity,tss)

	scaledAffinities={}
	if (args.signalScale != ""):
		scaledAffinities=extractTF_Affinity(genesInOpenChromatin,args.signalScale,tss,decay,loopsactivated, looptable)
		if (decay):
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Affinity_Gene_View.txt"),tss)
		else:
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Affinity_Gene_View.txt"),tss)	

main()
