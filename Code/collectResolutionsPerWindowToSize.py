import argparse
import datetime

import utils


def readIntraLoopsWithRes(loopsFile):
    lf = open(loopsFile, 'r')
    loops = {}

    for i in range(1, 23):
        loops[str(i)] = []

    loopID = 0
    for l in lf:
        s = l.split()
        if (len(s) >= 8):
            if (s[0] == s[3] and s[0] != 'X' and s[
                0] != 'Y'):  # check if loop is intra-chromosomal and not on X or Y chromsome
                loopID += 1
                resolution = int(s[2]) - int(s[1])
                loops[s[0]].append((loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7]), resolution))
    lf.close()
    return loops


def run(genes, intraLoops, windows):
    matchedLoops = {}

    for geneKey in genes:
        chrLoops = intraLoops[genes[geneKey][0]]
        pos = genes[geneKey][1]
        if not matchedLoops.has_key(geneKey):
            matchedLoops[geneKey] = {}

        for loop in chrLoops:
            resolution = loop[6]

            for window in windows:
                posRight = pos + window
                posLeft = pos - window

                if not matchedLoops[geneKey].has_key(window):
                    matchedLoops[geneKey][window] = {}

                if not matchedLoops[geneKey][window].has_key(resolution):
                    matchedLoops[geneKey][window][resolution] = []

                if (loop[1] >= posLeft and loop[1] <= posRight or loop[2] >= posLeft and loop[2] <= posRight):
                    matchedLoops[geneKey][window][resolution].append(loop)

                else:
                    if (loop[3] >= posLeft and loop[3] <= posRight or loop[4] >= posLeft and loop[4] <= posRight):
                        matchedLoops[geneKey][window][resolution].append(loop)

    return (matchedLoops)


def postProcessing(results, tss, intraLoops):
    header = 'window	count	resolution'
    output = ''

    for geneKey in results:
        windows = results[geneKey]

        for windowKey in windows:
            resolutions = windows[windowKey]

            for resolutionKey in resolutions:
                hits = resolutions[resolutionKey]
                output += str(windowKey) + '	'
                output += str(len(hits)) + '	'
                output += str(resolutionKey)
                output += '\n'

    now = datetime.datetime.now()
    filename = now.strftime("%Y-%m-%d-%H-%M") + '_ResolutionPerWindowToSize' + '.txt'

    utils.writeToFile(filename, header, output)


def collectResolutionsPerWindowToSize(annotationFile, loopsFile, windows):
    print 'Indexing TSS'
    tss = utils.readGTF(annotationFile)
    print 'Indexing Loops'
    intraLoops = readIntraLoopsWithRes(loopsFile)

    print 'Preprocessing'
    for window in windows:
        if (window < 100):
            print 'Window radius too small... please use greater values e.g. 100 and above.'
            return 1;

    print 'Running core algorithm'
    results = run(tss, intraLoops, windows)

    postProcessing(results, tss, intraLoops)

    return 0


######
##
##	Required arguments: 1) annotation File (.gtf), 2) loop file in Hi-C loop format, 3) a comma seperatated enumeration of window radii (window radii around TSS in thousand bases)
##
######

# preparation-routine 
parser = argparse.ArgumentParser(
    description='Collects all loops in a window around genes and extracts some statistical values')
parser.add_argument('annotation', help='Path to an annotation file')
parser.add_argument('loops', help='Path to a loop file')
parser.add_argument('windows', type=str, help='A list of window radii around the genestart which should be scanned')

args = parser.parse_args()
windows = args.windows.split(',')
windows = [int(numeric_string) for numeric_string in windows]

print 'Starting to collect data...'
collectResolutionsPerWindowToSize(args.annotation, args.loops, windows)
print '\n-> Completed all!'
