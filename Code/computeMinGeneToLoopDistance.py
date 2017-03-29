import argparse
import datetime

import utils


def run(tss, intraLoops, usemiddle, resolution):
    geneToLoopDistances = {}

    for geneKey in tss:
        chrLoops = intraLoops[tss[geneKey][0]]
        startpos = tss[geneKey][1]
        minDistance = -1
        loopres = -1
        for loop in chrLoops:
            if (loop[6] != resolution):
                continue

            if (usemiddle):
                middle = int(((float(loop[2]) - float(loop[1])) / 2) + float(loop[1]))
                if (minDistance == -1):
                    minDistance = abs(middle - startpos)
                    loopres = loop[6]

                if (abs(middle - startpos) < minDistance):
                    minDistance = abs(middle - startpos)
                    loopres = loop[6]

                middle = int(((float(loop[4]) - float(loop[3])) / 2) + float(loop[3]))
                if (abs(middle - startpos) < minDistance):
                    minDistance = abs(middle - startpos)
                    loopres = loop[6]

            else:
                if (minDistance == -1):
                    minDistance = abs(loop[1] - startpos)
                    loopres = loop[6]

                if (abs(loop[1] - startpos) < minDistance):
                    minDistance = abs(loop[1] - startpos)
                    loopres = loop[6]

                if (abs(loop[2] - startpos) < minDistance):
                    minDistance = abs(loop[2] - startpos)
                    loopres = loop[6]

                if (abs(loop[3] - startpos) < minDistance):
                    minDistance = abs(loop[3] - startpos)
                    loopres = loop[6]

                if (abs(loop[4] - startpos) < minDistance):
                    minDistance = abs(loop[4] - startpos)
                    loopres = loop[6]

        geneToLoopDistances[geneKey] = (minDistance, loopres)  # also append resolution
    return geneToLoopDistances


def postProcessing(geneToLoopDistances, usemiddle, resolution):
    header = 'geneID	distance	resolution'

    output = ''
    for geneKey in geneToLoopDistances:
        output += str(geneKey) + '	'
        output += str(geneToLoopDistances[geneKey][0]) + '	'
        output += str(geneToLoopDistances[geneKey][1])
        output += '\n'

    now = datetime.datetime.now()
    if (usemiddle):
        filename = now.strftime("%Y-%m-%d-%H-%M") + '_MinDistances_Res' + str(resolution) + '_Middle' + '.txt'
    else:
        filename = now.strftime("%Y-%m-%d-%H-%M") + '_MinDistances_Res' + str(resolution) + '.txt'

    utils.writeToFile(filename, header, output)
    return 0


# Main entry point for algorithm, requires:
#	an annotation-file
# 	a loop-file
def computeMinGeneToLoopDistance(annotationFile, loopsFile, resolution, usemiddle):
    print 'Indexing TSS'
    tss = utils.readGTF(annotationFile)
    print 'Indexing Loops'
    intraLoops = utils.readIntraLoops(loopsFile)

    print 'Running core algorithm'
    print 'Using middle: ' + str(usemiddle)
    results = run(tss, intraLoops, usemiddle, resolution)

    postProcessing(results, usemiddle, resolution)

    return 0


######
##
##	Required arguments: 1) annotation File (.gtf), 2) loop file in Hi-C loop format 3) Hi-C resolution in bp
##
######

# preparation-routine 
parser = argparse.ArgumentParser(description='Computes for each gene the distance to it\'s nearest loop.')
parser.add_argument('annotation', help='Path to an annotation file')
parser.add_argument('loops', help='Path to a loop file')
parser.add_argument('resolution', help='Only consider loops of the given resolution')
parser.add_argument('middle', help='Defines whether to use the middle or an edge of a loop for distance calculation.')

args = parser.parse_args()
usemiddle = False

if (args.middle.upper() == "TRUE"):
    usemiddle = True

print 'Starting to collect data...'
computeMinGeneToLoopDistance(args.annotation, args.loops, int(args.resolution), usemiddle)
print '\n-> Completed all!'
