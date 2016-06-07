import argparse
import datetime
import math
import utils


# Computing per gene TF affinities

# Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
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

# Reads the txt file containing TF-scores. Extracts the regions of open chromatin.
# They are returned as a dictionary(key: #chromosom, item:[(start,end)])
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


def extractTF_Affinity(openChromatinInGenes,filename,tss,expDecay,loopsactivated,loopOCregions,geneloops,looplookuptable,loopdecay):
	geneAffinities={}
	tfpa=readTFPA(filename)
	for geneID in tss:
		startpos=tss[geneID][1]
		if(openChromatinInGenes.has_key(geneID)):
			regions = openChromatinInGenes[geneID]
			for region in regions:
				if(tfpa.has_key(region)):
					s = tfpa[region]
					middles=s[0].split(":")[1].split("-")
					middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
					
					if (expDecay):
						factor=math.exp(-(float(float(abs(startpos-middle))/5000.0)))
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
		elif(loopsactivated):	# insert loop part here
			if(geneloops.has_key(geneID)):
				loops = geneloops[geneID]
				for tupel in loops:
					if(loopOCregions.has_key(tupel[0])):
						regions = loopOCregions[tupel[0]]
						aff = []
						
						for region in regions[0]:	# walk through left side of loop
							loci = tss[geneID][0]+":"+str(region[0])+"-"+str(region[1])
							s = tfpa[loci]
							middles=s[0].split(":")[1].split("-")
							if(tupel[1]):
								middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
							else:
								loopprops = looplookuptable[tupel[0]]
								distance = loopprops[3] - loopprops[2]
								middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))+distance
							aff.append((middle, s))
							
						for region in regions[1]:	#walk through right side of loop
							loci = tss[geneID][0]+":"+str(region[0])+"-"+str(region[1])
							s = tfpa[loci]
							middles=s[0].split(":")[1].split("-")
							if(not tupel[1]):
								middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))
							else:
								loopprops = looplookuptable[tupel[0]]
								distance = loopprops[3] - loopprops[2]
								middle=int(((float(middles[1])-float(middles[0]))/2)+float(middles[0]))-distance
							aff.append((middle, s))
						for afftupel in aff:
							if (loopdecay):
								factor=math.exp(-(float(float(abs(startpos-afftupel[0]))/5000.0)))
								if (geneAffinities.has_key(geneID)):
									geneAffinities[geneID]=aggregateAffinity(geneAffinities[geneID],afftupel[1][1:],factor)
								else:
									numbers=afftupel[1][1:]
									for i in xrange(0,len(numbers)-1):
										numbers[i]=float(factor)*float(numbers[i])
									geneAffinities[geneID]=numbers				
							else:
								if (geneAffinities.has_key(geneID)):
									geneAffinities[geneID]=aggregateAffinity(geneAffinities[geneID],afftupel[1][1:],1.0)
								else:
									geneAffinities[geneID]=afftupel[1][1:]
	
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
	
def createGenesWithLoopsFile(geneLoops, filename):
	header ="geneID"
	body = ''
	for geneID in geneLoops:
		if(len(geneLoops[geneID]) > 0):
			body += geneID
			body += '\n'
	utils.writeToFile(filename, header, body)
	
def intersectRegions(oc, loops):
	looptable = {}
	for chrKey in oc:
		for octupel in oc[chrKey]:
			if(loops.has_key(chrKey)):
				chrLoops = loops[chrKey]
				for loop in chrLoops:
					# open chromatin region in...
					#	[  ----  ]
					if(octupel[0] >= loop[1] and octupel[1] <= loop[2]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# --[------  ]
					elif(octupel[0] <= loop[1] and octupel[1] <= loop[2] and octupel[1] >= loop[1]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# 	[  ------]--	
					elif(octupel[0] >= loop[1] and octupel[0] <= loop[2] and octupel[1] >= loop[2]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# --[--------]--	
					elif(octupel[0] < loop[1] and octupel[1] > loop[2]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][0].append(octupel)
					# same procedure for second, counterpart loop-site
					#	[  ----  ]
					elif(octupel[0] >= loop[3] and octupel[1] <= loop[4]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# --[------  ]
					elif(octupel[0] <= loop[3] and octupel[1] <= loop[4] and octupel[1] >= loop[3]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# 	[  ------]--	
					elif(octupel[0] >= loop[3] and octupel[0] <= loop[4] and octupel[1] >= loop[4]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
					# --[--------]--	
					elif(octupel[0] < loop[3] and octupel[1] > loop[4]):
						if(not looptable.has_key(loop[0])):
							looptable[loop[0]] = ([], [])
						looptable[loop[0]][1].append(octupel)
						
	return looptable


def filterLoops(loops, resolution):
	for chrKey in loops:
		chrLoops = loops[chrKey]
		for loop in chrLoops:
			if((loop[2]-loop[1]) != resolution):
				chrLoops.remove(loop)
	
				
def findLoopsNearbyGenes(tss, loops, loopwindows, usemiddle):
	geneLoops = {}
	loopwindows = int(loopwindows/2)
	for geneID in tss:
		startpos = tss[geneID][1]
		if(loops.has_key(tss[geneID][0])):
			chrLoops = loops[tss[geneID][0]]
			if(not geneLoops.has_key(geneID)):
				geneLoops[geneID] = []
			for loop in chrLoops:
				foundleftone = False
				foundrightone = False
				if(usemiddle):
					
					leftmiddle = int(((float(loop[2])-float(loop[1]))/2)+float(loop[1]))
					if(abs(leftmiddle - startpos) <= loopwindows):
						foundleftone = True
						
					rightmiddle = int(((float(loop[3])-float(loop[4]))/2)+float(loop[3]))
					if (abs(rightmiddle - startpos) <= loopwindows):
						foundrightone = True
					
					if(foundleftone and foundrightone):
						if(abs(leftmiddle - startpos) <= abs(rightmiddle - startpos)):
							geneLoops[geneID].append((loop[0], True))
						else:
							geneLoops[geneID].append((loop[0], False))
						
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
						else:
							geneLoops[geneID].append((loop[0], False))
					
					elif(foundleftone):
						geneLoops[geneID].append((loop[0], True))
					
					elif(foundrightone):
						geneLoops[geneID].append((loop[0], False))
				
	return geneLoops
					

def main():
	parser=argparse.ArgumentParser(prog="annotateTSS.py")
	parser.add_argument("gtf",nargs=1,help="Genome annotation file")
	parser.add_argument("affinity",nargs=1,help="TRAP generated TF Affinity file")
	parser.add_argument("--geneViewAffinity",nargs="?",help="Name of the gene view affinity files. If this is not specified, the prefix of the input files will be used.",default="")
	parser.add_argument("--windows",nargs="?",help="Size of the considered window around the TSS. Default is 3000.",default=3000,type=int)
	parser.add_argument("--decay",nargs="?",help="True if exponential decay should be used, False otherwise. Default is True",default="True")
	parser.add_argument("--signalScale",nargs="?",help="If the name of the scaled affinity file is provied, a Gene view file is computed for those Affinity values.",default="")
	parser.add_argument("--loopfile",nargs="?",help="If the name of the loop file is provied, all open chromatin regions will be intersected with loop regions around the TSS of each gene.",default="")
	parser.add_argument("--loopwindows",nargs="?",help="Defines the window-size around the TSS in which all loops are considered for intersecting with openChromatin regions.",default=10000,type=int)
	parser.add_argument("--resolution",nargs="?",help="Defines the Hi-C resolution of the loops which should be considered. Uses the smallest one found if the a resolution is not available in the loopfile.",default=5000,type=int)
	parser.add_argument("--usemiddle",nargs="?",help="Defines whether to use the middle of a loop to decide if a loop lies withing a window or the edges.",default="False")
	parser.add_argument("--loopdecay",nargs="?",help="True if exponential decay should be used for oc regions in loops, False otherwise. Default is False",default="False")
	args=parser.parse_args() 

	prefixs=args.affinity[0].split(".")
	prefix=prefixs[0]
	if (args.geneViewAffinity==""):
		args.geneViewAffinity=prefix+"_Affinity_Gene_View.txt"

	if (args.decay.upper()=="FALSE") or (args.decay=="0"):
		decay=False
	else:
		decay=True
	
	if (args.loopdecay.upper()=="FALSE") or (args.decay=="0"):
		loopdecay=False
	else:
		loopdecay=True
		
	now = datetime.datetime.now()
	print 'Start time: ' + now.strftime("%Y-%m-%d-%H-%M-%S")
	
	# Check arguments
	
	# Extract TSS of GTF files
	tss=readGTF(args.gtf[0])
	# Load open chromatin positions from TF-Affinity file
	oC=readOC_Region(args.affinity[0])
	# Create a TF name index
	tfNames=tfIndex(args.affinity[0])
	shift=int(args.windows/2)
	
	loopwindows = args.loopwindows
	resolution = args.resolution
	loopsactivated = False
	loopOCregions = {}
	geneloops = {}
	looplookuptable = {}
	
	# Extract loops of Hi-C loopfile
	if(args.loopfile != ""):
		if(args.usemiddle.upper() == "TRUE"):
			usemiddle = True
		else:
			usemiddle = False
		loopsactivated = True
		
		loops = utils.readIntraLoops(args.loopfile)
		# filter loops and keep user-defined resolution only
		filterLoops(loops, resolution)
		geneloops = findLoopsNearbyGenes(tss, loops, loopwindows, usemiddle)
		looplookuptable = utils.readIntraLoopsLookupTable(args.loopfile) # maybe use this later on
		
		# intersect openchromatin regions with all loopregions in TSS windows
		intersectResults = intersectRegions(oC, loops)
		#oC = intersectResults[0]
		loopOCregions = intersectResults
		
	
	# Determine gene windows in open chromatin regions
	openChromatinInGenes={}
	for gene in tss.keys():
		if (oC.has_key(tss[gene][0])):
			for tupel in oC[tss[gene][0]]:
				
				if(loopsactivated):
					found = False
					if(geneloops.has_key(gene)):
						for geneloop in geneloops[gene]:
							if(loopOCregions.has_key(geneloop[0])):
								regions = loopOCregions[geneloop[0]]
								if(geneloop[1] and regions[0].count(tupel) > 0):
									found = True
									break
								elif((not geneloop[1]) and regions[1].count(tupel) > 0):
									found = True
									break
					if(found):
						continue
				
				loci = tss[gene][0]+":"+str(tupel[0])+"-"+str(tupel[1])
				# Right border of window <= Right border of open chromatin
				if (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift >= tupel[0]):
					# Left border of window >= Left border of open chromatin ==> Window inside open chromatin
					if (openChromatinInGenes.has_key(gene)):		
						openChromatinInGenes[gene]+=[loci]
					else:
						openChromatinInGenes[gene]=[loci]					
				# Right border of window >= Left border of open chromatin ==> Window enters open chromatin in the 3' end and stops in the tss window
				elif (tss[gene][1]+shift <= tupel[1]) and (tss[gene][1]-shift < tupel[0]) and (tss[gene][1]+shift > tupel[0]):
					if (openChromatinInGenes.has_key(gene)):		
						openChromatinInGenes[gene]+=[loci]
					else:
						openChromatinInGenes[gene]=[loci]					
				# Right border of window > Right border of open chromatin
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift < tupel[0]):
					# Left border of window <= Left border of open chromatin ==> Window is larger than open chromatin
					if (openChromatinInGenes.has_key(gene)):		
						openChromatinInGenes[gene]+=[loci]
					else:
						openChromatinInGenes[gene]=[loci]
				# Left border of window <= Right border of open chromain ==> Window enters open chromatin in the 5' end stops in the tss window
				elif (tss[gene][1]+shift > tupel[1]) and (tss[gene][1]-shift >= tupel[0]) and (tss[gene][1]-shift < tupel[1]):
					if (openChromatinInGenes.has_key(gene)):		
						openChromatinInGenes[gene]+=[loci]
					else:
						openChromatinInGenes[gene]=[loci]
	
	# Extract bound transcription factors
	affinities=extractTF_Affinity(openChromatinInGenes,args.affinity[0],tss,decay,loopsactivated,loopOCregions,geneloops,looplookuptable,loopdecay)
	
	if (decay):
		createAffinityFile(affinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Affinity_Gene_View.txt"),tss)	
	else:
		createAffinityFile(affinities,tfNames,args.geneViewAffinity,tss)

	scaledAffinities={}
	if (args.signalScale != ""):
		scaledAffinities=extractTF_Affinity(openChromatinInGenes,args.signalScale,tss,decay,loopsactivated,loopOCregions,geneloops,looplookuptable,loopdecay)
		if (decay):
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Decay_Scaled_Affinity_Gene_View.txt"),tss)
		else:
			createAffinityFile(scaledAffinities,tfNames,args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_Scaled_Affinity_Gene_View.txt"),tss)
	
	if(loopsactivated):
		createGenesWithLoopsFile(geneloops, args.geneViewAffinity.replace("_Affinity_Gene_View.txt","_GenesWithLoops.txt"))
	
	now = datetime.datetime.now()
	print 'End time: ' + now.strftime("%Y-%m-%d-%H-%M-%S")

main()
