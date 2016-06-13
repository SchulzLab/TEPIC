import argparse
import utils

def filterInclusions(intraLoops, inclusions):
	for chrKey in inclusions:
		chrInclusions = inclusions[chrKey]
		for inclusion in chrInclusions:
			if(inclusion[0][6] < inclusion[1][6]):
				intraLoops[chrKey].remove(inclusion[1])
			elif(inclusion[0][6] > inclusion[1][6]):
				intraLoops[chrKey].remove(inclusion[0])
	return

def run(intraLoops):
	inclusions = {}
	for chrKey in intraLoops:
		if(not inclusions.has_key(chrKey)):
			inclusions[chrKey] = []
		chrLoops = intraLoops[chrKey]
		for loop in chrLoops:
			for otherloop in chrLoops:
				if(loop[0] == otherloop[0]):# same loop
					continue
				if(loop[6] <= otherloop[6]):# same resolution or lower
					continue
				if(loop[1] <= otherloop[1] and otherloop[2] <= loop[2] and loop[3] <= otherloop[3] and otherloop[4] <= loop[4]):
					inclusions[chrKey].append((loop, otherloop))
	return inclusions


def postProcessing(results):
	header = ''
	output = ''

	for chrKey in results:
		chrInclusions = results[chrKey]
		if(len(chrInclusions) == 0):
			output += 'No inclusions on Chr.' + str(chrKey) + ' found'
			output += '\n'
		else:
			output += 'Inclusions on Chr.' + str(chrKey) + ':'
			output += '\n'
			for looppair in chrInclusions:
				output += str(looppair[0])
				output += str(looppair[1])
				output += '\n'

	utils.writeToFile('detectionoutput.txt', header, output)
	return 0


# Main entry point for algorithm, requires:
# 	a loop-file
def detectLoopInclusions(loopsFile):
	print 'Indexing Loops'
	intraLoops = utils.readIntraLoops(loopsFile)

	print 'Running core algorithm'
	results = run(intraLoops)

	postProcessing(results)

	return 0


######
##
##	Required arguments: 1) loop file in Hi-C loop format
##
######
def main():
	# preparation-routine
	parser = argparse.ArgumentParser(description='Checks whether there are (Intrachromosomal) loops of different resolutions in a loop file, so that a loop of a lower resolution contains a loop with a higher resolution.')
	parser.add_argument('loops', help='Path to a loop file')

	args = parser.parse_args()


	print 'Starting to collect data...'
	detectLoopInclusions(args.loops)
	print '\n-> Completed all!'
