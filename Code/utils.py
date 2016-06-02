

# Reads an annotation file and returns a dictionary:
# {key : value} -> {gene : (chr., starting position)}
def readGTF(annotationFile):
	gtf=open(annotationFile, "r")
	tss={}
	for l in gtf:
		s=l.split()
		if (len(s) >=9):
			if (s[2]=="gene" and s[0]!='chrM' and s[0]!='chrX' and s[0]!='chrY'):	
				if (s[6]=="+"):
					tss[s[9].replace(";","")]=(s[0].replace("chr",""),int(s[3]))
				else:
					tss[s[9].replace(";","")]=(s[0].replace("chr",""),int(s[4]))
	gtf.close()
	return tss


# Reads a loop file and returns a dictionary containing all intrachromosomal loops per chromosome in a list:
# {key : value} -> {chr : [(loopID, start X, end X, start Y, end Y, observations count, resolution)]}
def readIntraLoops(loopsFile):
	lf = open(loopsFile, 'r')
	loops = {}
	
	loopID = 0
	for l in lf:
		s=l.split()
		if (len(s) >=8):
			if(s[0] == s[3]): # check if loop is intra-chromosomal
				loopID += 1
				resolution = int(s[2]) - int(s[1])
				if(not loops.has_key(s[0])):
					loops[s[0]] = []
				loops[s[0]].append((loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7]), resolution))
	lf.close()
	return loops

# Reads a loop file and returns a dictionary containg all intrachromosomal loops accessable by an unique ID:
# {key : value} -> {loopID : (loopID, start X, end X, start Y, end Y, observations count, resolution)}
# we keep the loop format for consistent programming and accessing tupel fields
def readIntraLoopsLookupTable(loopsFile):
	lf = open(loopsFile, 'r')
	loops = {}
	
	loopID = 0
	for l in lf:
		s=l.split()
		if (len(s) >=8):
			if(s[0] == s[3] and s[0]!='X' and s[0]!='Y'): # check if loop is intra-chromosomal
				loopID += 1
				resolution = int(s[2]) - int(s[1])
				loops[loopID] = (loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7]), resolution)
	lf.close()
	return loops


# Reads a loop file and returns a dictionary containg all intrachromosomal loops per chromosome in a list:
#TODO: define format	
def readInterLoops(loopsFile):
	#TODO: implement this
	print 'Not implemented yet!'
	return 0


# Writes given header and body to a file defined by filename
def writeToFile(filename, header, body ):
	print 'Writing contents of ' + filename + ' to disk...' 
	f = open(filename, 'w')
	f.write(str(header))
	f.write('\n')
	f.write(str(body))
	f.close()
	print 'Finished'
