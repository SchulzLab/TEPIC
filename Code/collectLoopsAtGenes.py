import sys
import datetime
import argparse


# Reads an annotation file and returns a dictionary:
# {key : value} -> {gene : (chr., starting position)}
def readGTF(annotationFile):
	gtf=open(annotationFile, "r")
	tss={}
	for l in gtf:
		s=l.split()
		if (len(s) >=9):
			if (s[2]=="gene" and s[0]!='chrM'):	
				if (s[6]=="+"):
					tss[s[9]]=(s[0].replace("chr",""),int(s[3]))
				else:
					tss[s[9]]=(s[0].replace("chr",""),int(s[4]))
	gtf.close()
	return tss


# Reads a loop file and returns a dictionary containg all intrachromosomal loops per chromosome in a list:
# {key : value} -> {chr : [(loopID, start X, end X, start Y, end Y, observations count)]}
def readIntraLoops(loopsFile):
	lf = open(loopsFile, 'r')
	loops = {}
	
	for i in range(1, 23):
		loops[str(i)] = []
	loops['X'] = []
	loops['Y'] = []
	
	loopID = 0
	for l in lf:
		s=l.split()
		if (len(s) >=8):
			if(s[0] == s[3]): # check if loop is intra-chromosomal
				loopID += 1
				loops[s[0]].append((loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7])))
	lf.close()
	return loops


# Reads a loop file and returns a dictionary containg all intrachromosomal loops per chromosome in a list:
#TODO: define format	
def readInterLoops(loopsFile):
	#TODO: implement this
	print 'Not implemented yet!'
	return 0


# Iterates over all genes and for each gene over all loops of its respective chromsome to find loops within a given radius of the TSS
# Returns:
# 	all matched loops
#		{key : value} -> {gene : [loopID]}
# 	the loopset
#		set of { loopID }
# 	all loops per chromosome
#		{key : value} -> {chr : set of {(loopID, start X, end X, start Y, end Y, observations count)}}
def run(genes, loops, radius):
	matchedLoops = {}
	loopSet = set()
	loopsPerChromosome = {}
	
	for geneKey in genes:
		chrLoops = loops[genes[geneKey][0]]
		pos = genes[geneKey][1]
		posRight = pos + radius
		posLeft = pos -radius
		
		if not matchedLoops.has_key(geneKey):
			matchedLoops[geneKey] = []
			
		if not loopsPerChromosome.has_key(genes[geneKey][0]):
			loopsPerChromosome[genes[geneKey][0]] = set()
			
		for loop in chrLoops:
			
			if(loop[1] >= posLeft and loop[1] <= posRight or loop[2] >= posLeft and loop[2] <= posRight):
				
				matchedLoops[geneKey].append(loop[0])
				loopSet.add(loop)
				loopsPerChromosome[genes[geneKey][0]].add(loop)
				 
			else:
				if(loop[3] >= posLeft and loop[3] <= posRight or loop[4] >= posLeft and loop[4] <= posRight):
					
					matchedLoops[geneKey].append(loop[0])
					loopSet.add(loop)
					loopsPerChromosome[genes[geneKey][0]].add(loop)
			
		
	return (matchedLoops, loopSet, loopsPerChromosome)


# Writes given header and body to a file defined by filename
def writeToFile(filename, header, body ):
	print 'Writing contents of ' + filename + ' to disk...' 
	f = open(filename, 'w')
	f.write(str(header))
	f.write('\n')
	f.write(str(body))
	f.close()
	print 'Finished'

# Filters the 'X' top genes which have the most loops around their TSS	
def filterTopXGenes(matchedLoops, loopSet, maxl):
	#TODO: implement this
	print 'not implemented'
	return 0

# Aggregates various statistical properties from the results and writes it to file
def postProcessing(results, tss, intraLoops, radius):
	
	totalGeneCount = len(tss)
	
	totalLoopCount = 0
	for chrLoops in intraLoops:
		totalLoopCount += len(intraLoops[chrLoops])
	
	matchedLoops = results[0]
	loopSet = results[1]
	loopsPerChromosome = results[2]
	
	hitCounter = 0
	maxl = 0
	
	output = {}
	header = 'ChrNo	Counts	Remaining	Total	Coverage'
	outstring = ''
	
	print 'Extracting statistics and preparing output'
	for geneKey in matchedLoops:
		hits = len(matchedLoops[geneKey])
		if hits > maxl:
			maxl = hits
		hitCounter += hits
	
	for chrKey in loopsPerChromosome:
		outputKey = chrKey
		if not output.has_key(outputKey):
			output[outputKey] = 0
		output[outputKey] += len(loopsPerChromosome[chrKey])
	
	for chrKey in output:
		outstring += str(chrKey) + '	'
		outstring += str(output[chrKey]) + '	'
		outstring += str(len(intraLoops[chrKey]) - output[chrKey]) + '	'
		outstring += str(len(intraLoops[chrKey])) + '	'
		if(len(intraLoops[chrKey]) > 0):
			outstring += str('%.1f' % round(float(len(loopsPerChromosome[chrKey]))/float(len(intraLoops[chrKey])) * 100, 1))
		else:
			outstring += str(0.0)
		outstring += '\n'
	
	outstring += 'ALL	'
	outstring += str(sum(output.values())) + '	'
	outstring += str(sum(len(v) for v in intraLoops.itervalues()) - len(loopSet)) + '	'
	outstring += str(sum(len(v) for v in intraLoops.itervalues())) + '	'
	outstring += str('%.1f' % round(float(len(loopSet))/float(totalLoopCount) * 100, 1))
	outstring += '\n'
	
	print '\n###############'
	print 'Parsed ' + str(totalGeneCount) + ' genes in total.' 
	print 'Parsed ' + str(totalLoopCount) + ' intrachromosomal loops in total.'
	print 'Found ' + str(len(loopSet)) + ' loops within a distance of ' + str(radius) + ' bases around the annotated TSS.'
	print 'That are ' +  str('%.1f' % round(float(len(loopSet))/float(totalLoopCount) * 100, 1)) + '% of all loops.'
	print 'Observed a maximum of ' + str(maxl) + ' loops nearby a single gene.' 
	print '###############\n'
	
	now = datetime.datetime.now()
	filename = now.strftime("%Y-%m-%d-%H-%M") + '_' + str(radius/1000) + 'kb' + '.txt'
	
	writeToFile(filename, header, outstring)
	
	return 0
		

# Main entry point for algorithm, requires:
#	an annotation-file
# 	a loop-file
# 	a radius in 1000 bases, which is >= 100	
def collectLoopsAtGenes(annotationFile, loopsFile, radius):
	print 'Indexing TSS'
	tss = readGTF(annotationFile)
	print 'Indexing Loops'
	intraLoops = readIntraLoops(loopsFile)
	
	if (radius < 100):
		print 'Radius too small... please use greater values e.g. 100 and above.'
		return 1;
	
	print 'Running core algorithm'
	results = run(tss, intraLoops, radius)
	
	postProcessing(results, tss, intraLoops, radius)
	
	return 0


######
##
##	Required arguments: 1) annotation File (.gtf), 2) loop file in Hi-C loop format, 3) radius (window radius around TSS in thousand bases)
##
######

# preparation-routine 
parser = argparse.ArgumentParser(description='Collects all loops in a window around genes and extracts some statistical values')
parser.add_argument('annotation', help='Path to an annotation file')
parser.add_argument('loops', help='Path to a loop file')
parser.add_argument('radius', type=int, help='A radius around the genestart which should be scanned')

args=parser.parse_args()

print 'Starting to collect data...'
collectLoopsAtGenes(args.annotation, args.loops, args.radius)
print '\n-> Completed all!'
