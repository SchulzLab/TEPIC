import argparse
import datetime
import math

import detectLoopInclusions
import utils


# Computing per gene TF affinities


# Reads the txt file containing TF-scores. Extracts the regions of open chromatin.
# They are returned as a dictionary(key: #chromosom, item:[(start,end)])
def readOC_Region(filename):
    tfpa = open(filename, "r")
    oC = {}
    for l in tfpa:
        s = l.split()[0]
        ds = s.split(":")
        if len(ds) >= 2:
            se = ds[1].split("-")
            if ds[0].replace("chr", "") not in oC:
                oC[ds[0].replace("chr", "")] = [(int(se[0]), int(se[1]))]
            else:
                oC[ds[0].replace("chr", "")] += [(int(se[0]), int(se[1]))]
    tfpa.close()
    return oC


def readTFPA(filename):
    f = open(filename, "r")
    f.readline()
    tfpa = {}
    for l in f:
        s = l.split()
        tfpa[s[0]] = s
    f.close()
    return tfpa


def aggregateAffinity(old, new, factor):
    for i in xrange(0, len(old)):
        old[i] = float(old[i]) + factor * float(new[i])
    return old


def extractTF_Affinity(openChromatinInGenes, filename, tss, expDecay, loopsactivated, loopOCregions, geneloops,
                       looplookuptable, loopdecay, loopcountscaling, countersiteonly, doublefeatures, enhancer):
    geneAffinities = {}
    loopAffinities = {}
    tfpa = readTFPA(filename)
    for geneID in tss:
        startpos = tss[geneID][1]
        if geneID in openChromatinInGenes:
            regions = openChromatinInGenes[geneID]
            for region in regions:
                if region in tfpa:
                    s = tfpa[region]
                    middles = s[0].split(":")[1].split("-")
                    middle = int(((float(middles[1]) - float(middles[0])) / 2) + float(middles[0]))

                    if expDecay:
                        factor = math.exp(-(float(float(abs(startpos - middle)) / 5000.0)))
                        if geneID in geneAffinities:
                            geneAffinities[geneID] = aggregateAffinity(geneAffinities[geneID], s[1:], factor)
                        else:
                            numbers = s[1:]
                            for i in xrange(0, len(numbers)):
                                numbers[i] = float(factor) * float(numbers[i])
                            geneAffinities[geneID] = numbers
                    else:
                        if geneID in geneAffinities:
                            geneAffinities[geneID] = aggregateAffinity(geneAffinities[geneID], s[1:], 1.0)
                        else:
                            geneAffinities[geneID] = s[1:]
        if loopsactivated:
            chrKey = tss[geneID][0]
            if geneID in geneloops:
                loops = geneloops[geneID]
                for tupel in loops:
                    if tupel[0] in loopOCregions:
                        left_side_closer = tupel[1]
                        regions = loopOCregions[tupel[0]]
                        aff = []
                        for region in regions[0]:  # walk through left side of loop
                            if countersiteonly and left_side_closer:
                                break
                            if not left_side_closer and not OCoverlapswithEnhancer(chrKey, region, enhancer):
                                continue
                            loci = tss[geneID][0] + ":" + str(region[0]) + "-" + str(region[1])
                            if loci in tfpa:
                                s = tfpa[loci]
                                middles = s[0].split(":")[1].split("-")
                                middle = int(((float(middles[1]) - float(middles[0])) / 2) + float(middles[0]))
                                loopprops = looplookuptable[tupel[0]]
                                loopcount = loopprops[5]  # extract observation count of loop
                                # check if current loop-site is further away from the gene than other loop-site
                                if not tupel[1]:
                                    distance = loopprops[3] - loopprops[2]
                                    middle += distance
                                aff.append((middle, s, loopcount))

                        for region in regions[1]:  # walk through right side of loop
                            if countersiteonly and not left_side_closer:
                                break
                            if left_side_closer and not OCoverlapswithEnhancer(chrKey, region, enhancer):
                                continue
                            loci = tss[geneID][0] + ":" + str(region[0]) + "-" + str(region[1])
                            if loci in tfpa:
                                s = tfpa[loci]
                                middles = s[0].split(":")[1].split("-")
                                middle = int(((float(middles[1]) - float(middles[0])) / 2) + float(middles[0]))
                                loopprops = looplookuptable[tupel[0]]
                                loopcount = loopprops[5]  # extract observation count of loop
                                # check if current loop-site is further away from the gene than other loop-site
                                if tupel[1]:
                                    distance = loopprops[3] - loopprops[2]
                                    middle -= distance
                                aff.append((middle, s, loopcount))

                        if not doublefeatures:
                            loopAffinities = geneAffinities
                        for afftupel in aff:
                            factor = 1.0
                            if loopdecay:
                                factor = math.exp(-(float(float(abs(startpos - afftupel[0])) / 5000.0)))
                            if loopcountscaling:
                                factor = factor * afftupel[2]
                            if geneID in loopAffinities:
                                loopAffinities[geneID] = aggregateAffinity(loopAffinities[geneID], afftupel[1][1:],
                                                                           factor)
                            else:
                                numbers = afftupel[1][1:]
                                for i in xrange(0, len(numbers)):
                                    numbers[i] = float(factor) * float(numbers[i])
                                loopAffinities[geneID] = numbers

    if loopsactivated and not doublefeatures:
        return geneAffinities, {}

    return geneAffinities, loopAffinities


def tfIndex(filename):
    tfpa = open(filename, "r")
    l = tfpa.readline()
    tfpa.close()
    return l.split()


def createAffinityFile(affs, tfNames, filename, tss, loopsactivated, doublefeatures):
    output = open(filename, "w")
    header = "geneID"
    for element in tfNames:
        header += '\t' + str(element)
    if loopsactivated and doublefeatures:
        for element in tfNames:
            header += '\t' + str(element) + '_LOOP'
    output.write(header + '\n')
    affinities = affs[0]
    loopaffinities = affs[1]
    for Gene in tss.keys():
        line = str(Gene.replace("\"", "").replace(";", "").split(".")[0])
        if Gene in affinities:
            for entry in affinities[Gene]:
                line += '\t' + str(entry)
        else:
            for i in xrange(0, len(tfNames)):
                line += '\t' + str(0.0)
        if loopsactivated and doublefeatures:
            if Gene in loopaffinities:
                for entry in loopaffinities[Gene]:
                    line += '\t' + str(entry)
            else:
                for i in xrange(0, len(tfNames)):
                    line += '\t' + str(0.0)
        output.write(line + '\n')
    output.close()


def createSparseFile(affinities, tfNames, filename, tss, number):
    if len(tfNames) < number:
        number = len(tfNames)
        print("Warning: The value of sparseRep is to large, representation will contain all possible TFs")
    output = open(filename, "w")
    header = "geneID\tTF\tAffinity\n"
    output.write(header)
    for Gene in tss.keys():
        tfList = []
        if Gene in affinities:
            geneID = str(Gene.replace("\"", "").replace(";", "").split(".")[0])
            temp = affinities[Gene]
            for i in xrange(0, len(tfNames)):
                tfList += [(tfNames[i], float(temp[i]))]
            tfList.sort(key=lambda tup: tup[1], reverse=True)
            for i in xrange(0, number):
                output.write(str(geneID) + "\t" + str(tfList[i][0]) + "\t" + str(tfList[i][1]) + "\n")
    output.close()


def makeTupels(values, names):
    l = []
    for i in xrange(0, len(values)):
        l += [(names[i], values[i])]
    return l


def createGenesWithLoopsFile(geneLoops, filename):
    header = "geneID"
    body = ''
    for geneID in geneLoops:
        if len(geneLoops[geneID]) > 0:
            body += geneID
            body += '\n'
    utils.writeToFile(filename, header, body)


# Returns True if the given OC regions lies within an annotated enhancer region
# In case of missing enhancer loadings, the methods returns True anyway.
def OCoverlapswithEnhancer(chrKey, octupel, enhancer):
    if not enhancer:
        return True

    if chrKey in enhancer:
        for e in enhancer[chrKey]:
            # open chromatin region in...
            # [  ----  ]
            if octupel[0] >= e[0] and octupel[1] <= e[1]:
                return True
            # --[------  ]
            elif octupel[0] <= e[0] <= octupel[1] <= e[1]:
                return True
            # [  ------]--
            elif e[0] <= octupel[0] <= e[1] <= octupel[1]:
                return True
            # --[--------]--
            elif octupel[0] < e[0] and octupel[1] > e[1]:
                return True

    return False


def intersectRegions(oc, loops):
    looptable = {}
    for chrKey in oc:
        for octupel in oc[chrKey]:
            if chrKey in loops:
                chrLoops = loops[chrKey]
                for loop in chrLoops:
                    # open chromatin region in...
                    # [  ----  ]
                    if octupel[0] >= loop[1] and octupel[1] <= loop[2]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # --[------  ]
                    elif octupel[0] <= loop[1] <= octupel[1] <= loop[2]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # [  ------]--
                    elif loop[1] <= octupel[0] <= loop[2] <= octupel[1]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # --[--------]--
                    elif octupel[0] < loop[1] and octupel[1] > loop[2]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # same procedure for second, counterpart loop-site
                    # [  ----  ]
                    elif octupel[0] >= loop[3] and octupel[1] <= loop[4]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # --[------  ]
                    elif octupel[0] <= loop[3] <= octupel[1] <= loop[4]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # [  ------]--
                    elif loop[3] <= octupel[0] <= loop[4] <= octupel[1]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # --[--------]--
                    elif octupel[0] < loop[3] and octupel[1] > loop[4]:
                        if loop[0] not in looptable:
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)

    return looptable


def filterLoops(loops, resolution):
    for chrKey in loops:
        chrLoops = loops[chrKey]
        for loop in chrLoops:
            if (loop[2] - loop[1]) != resolution:
                chrLoops.remove(loop)


# Find all loops belonging to a gene that are in the range of the given loopwindow around the TSS
# Return dictionary: key: geneID, value: [(loopID, BOOLEAN)], the boolean value is True if the left side
# of the loop is closer to the annotated tss, false if the right side is closer.
def findLoopsNearbyGenes(tss, loops, loopwindows, usemiddle):
    geneLoops = {}
    loopwindows = int(loopwindows / 2)
    for geneID in tss:
        startpos = tss[geneID][1]
        if tss[geneID][0] in loops:
            chrLoops = loops[tss[geneID][0]]
            if geneID not in geneLoops:
                geneLoops[geneID] = []
            for loop in chrLoops:
                foundleftone = False
                foundrightone = False
                if usemiddle:

                    leftmiddle = int(((float(loop[2]) - float(loop[1])) / 2) + float(loop[1]))
                    if abs(leftmiddle - startpos) <= loopwindows:
                        foundleftone = True

                    rightmiddle = int(((float(loop[3]) - float(loop[4])) / 2) + float(loop[3]))
                    if abs(rightmiddle - startpos) <= loopwindows:
                        foundrightone = True

                    if foundleftone and foundrightone:
                        if abs(leftmiddle - startpos) <= abs(rightmiddle - startpos):
                            geneLoops[geneID].append((loop[0], True))
                        else:
                            geneLoops[geneID].append((loop[0], False))

                    elif foundleftone:
                        geneLoops[geneID].append((loop[0], True))

                    elif foundrightone:
                        geneLoops[geneID].append((loop[0], False))

                else:
                    leftpos = 0
                    if (abs(loop[1] - startpos) <= loopwindows) or (abs(loop[2] - startpos) <= loopwindows):
                        foundleftone = True
                        if abs(loop[1] - startpos) <= abs(loop[2] - startpos):
                            leftpos = abs(loop[1] - startpos)
                        else:
                            leftpos = abs(loop[2] - startpos)

                    rightpos = 0
                    if (abs(loop[3] - startpos) <= loopwindows) or (abs(loop[4] - startpos) <= loopwindows):
                        foundrightone = True
                        if abs(loop[3] - startpos) <= abs(loop[4] - startpos):
                            rightpos = abs(loop[3] - startpos)
                        else:
                            rightpos = abs(loop[4] - startpos)

                    if foundleftone and foundrightone:
                        if leftpos <= rightpos:
                            geneLoops[geneID].append((loop[0], True))
                        else:
                            geneLoops[geneID].append((loop[0], False))

                    elif foundleftone:
                        geneLoops[geneID].append((loop[0], True))

                    elif foundrightone:
                        geneLoops[geneID].append((loop[0], False))

    return geneLoops


def main():
    parser = argparse.ArgumentParser(prog="annotateTSS.py")
    parser.add_argument("gtf", nargs=1, help="Genome annotation file")
    parser.add_argument("affinity", nargs=1, help="TRAP generated TF Affinity file")
    # noinspection PyPep8
    parser.add_argument("--geneViewAffinity", nargs="?",
                        help="Name of the gene view affinity files. If this is not specified, the prefix of the input files will be used.",
                        default="")
    parser.add_argument("--windows", nargs="?", help="Size of the considered window around the TSS. Default is 3000.",
                        default=3000, type=int)
    parser.add_argument("--decay", nargs="?",
                        help="True if exponential decay should be used, false otherwise. Default is True.",
                        default="True")
    parser.add_argument("--signalScale", nargs="?",
                        help="If the name of a scaled affinity file is provided, a second GeneView file is computed.",
                        default="")
    # noinspection PyPep8
    parser.add_argument("--loopfile", nargs="?",
                        help="If the name of the Hi-C loop file is provided, all open chromatin regions will be intersected with loop regions around the TSS of each gene. Activates the Hi-C mode.",
                        default="")
    # noinspection PyPep8
    parser.add_argument("--loopwindows", nargs="?",
                        help="Defines the window-size around the TSS in which all loops are considered for intersecting with openChromatin regions.",
                        default=10000, type=int)
    parser.add_argument("--resolution", nargs="?",
                        help="Defines the Hi-C resolution of the loops which should be considered.", default=5000)
    # noinspection PyPep8
    parser.add_argument("--usemiddle", nargs="?",
                        help="Defines whether to use the middle of a loop to decide if a loop lies within a window or the edges.",
                        default="False")
    # noinspection PyPep8
    parser.add_argument("--loopdecay", nargs="?",
                        help="Set True if exponential decay should be used for oc regions in loops, false otherwise. Default is False.",
                        default="False")
    # noinspection PyPep8
    parser.add_argument("--loopcountscaling", nargs="?",
                        help="Set True if open chromatin regions inside loop-sites should be scaled with the loopcount given by the Hi-C loop-file, false otherwise. Default is False.",
                        default="False")
    # noinspection PyPep8
    parser.add_argument("--countersiteonly", nargs="?",
                        help="Set True if only the counter loop-site should affect the scores. The loop-site lying near the TSS will be disabled. Default is false.",
                        default="False")
    # noinspection PyPep8
    parser.add_argument("--doublefeatures", nargs="?",
                        help="Set True if the transcription factor feature space should be doubled in Hi-C mode. The Gene-TF-score will then be generated twice, once for the TSS associated TF score and once for the Hi-C loops. Default is false.",
                        default="False")
    parser.add_argument("--sparseRep", nargs="?",
                        help="Number of top TFs that should be contained in the sparse representation", default=0,
                        type=int)
    # noinspection PyPep8
    parser.add_argument("--enhancerfile", nargs="?",
                        help="Intersect Open Chromatin regions inside Hi-C loop regions with enhancer annotations. Only available in Hi-C mode.",
                        default="")
    args = parser.parse_args()

    # Check arguments
    prefixs = args.affinity[0].split(".")
    prefix = prefixs[0]
    if args.geneViewAffinity == "":
        args.geneViewAffinity = prefix + "_Affinity_Gene_View.txt"

    decay = False
    if (args.decay.upper() == "TRUE") or (args.decay == "1"):
        decay = True

    loopdecay = False
    if (args.loopdecay.upper() == "TRUE") or (args.loopdecay == "1"):
        loopdecay = True

    loopcountscaling = False
    if (args.loopcountscaling.upper() == "TRUE") or (args.loopcountscaling == "1"):
        loopcountscaling = True

    countersiteonly = False
    if (args.countersiteonly.upper() == "TRUE") or (args.countersiteonly == "1"):
        countersiteonly = True

    usemiddle = False
    if (args.usemiddle.upper() == "TRUE") or (args.usemiddle == "1"):
        usemiddle = True

    doublefeatures = False
    if (args.doublefeatures.upper() == "TRUE") or (args.doublefeatures == "1"):
        doublefeatures = True

    now = datetime.datetime.now()
    print 'Start time: ' + now.strftime("%Y-%m-%d-%H-%M-%S")

    # Extract TSS of GTF files
    tss = utils.readGTF(args.gtf[0])
    # Load open chromatin positions from TF-Affinity file
    oC = readOC_Region(args.affinity[0])
    # Create a TF name index
    tfNames = tfIndex(args.affinity[0])
    shift = int(args.windows / 2)
    # Initialize variables
    loopwindows = args.loopwindows
    resolution = args.resolution
    loopsactivated = False
    loopOCregions = {}
    geneloops = {}
    looplookuptable = {}
    enhancer = {}

    # Extract loops of Hi-C loopfile
    if args.loopfile != "":
        print 'Running annotation in Hi-C mode'
        loopsactivated = True

        loops = utils.readIntraLoops(args.loopfile)

        if str(resolution).upper() == "ALL":
            # detect inclusions of loops between different resolutions and remove them
            # keep all other loops of any resolution
            inclusions = detectLoopInclusions.run(loops)
            detectLoopInclusions.filterInclusions(loops, inclusions)
        else:
            # filter loops and keep user-defined resolution only
            resolution = int(args.resolution)
            filterLoops(loops, resolution)

        geneloops = findLoopsNearbyGenes(tss, loops, loopwindows, usemiddle)
        looplookuptable = utils.readIntraLoopsLookupTable(args.loopfile)

        # Read enhancer annotations
        if args.enhancerfile != "":
            enhancer = utils.readEnhancer(args.enhancerfile)

        # intersect openchromatin regions with all loopregions in TSS windows
        loopOCregions = intersectRegions(oC, loops)

    else:
        print 'Running annotation in original mode'

    # Determine gene windows in open chromatin regions
    openChromatinInGenes = {}
    for gene in tss.keys():
        if tss[gene][0] in oC:
            for tupel in oC[tss[gene][0]]:

                if loopsactivated:
                    found = False
                    if gene in geneloops:
                        for geneloop in geneloops[gene]:
                            if geneloop[0] in loopOCregions:
                                regions = loopOCregions[geneloop[0]]
                                if geneloop[1] and regions[0].count(tupel) > 0:
                                    found = True
                                    break
                                elif (not geneloop[1]) and regions[1].count(tupel) > 0:
                                    found = True
                                    break
                    if found:
                        continue

                loci = tss[gene][0] + ":" + str(tupel[0]) + "-" + str(tupel[1])
                # Right border of window <= Right border of open chromatin
                if (tss[gene][1] + shift <= tupel[1]) and (tss[gene][1] - shift >= tupel[0]):
                    # Left border of window >= Left border of open chromatin ==> Window inside open chromatin
                    if gene in openChromatinInGenes:
                        openChromatinInGenes[gene] += [loci]
                    else:
                        openChromatinInGenes[gene] = [loci]
                # Right border of window >= Left border of open chromatin ==> Window enters open chromatin in the 3'
                # end and stops in the tss window
                elif (tss[gene][1] + shift <= tupel[1]) and (tss[gene][1] - shift < tupel[0]) and (
                                tss[gene][1] + shift > tupel[0]):
                    if gene in openChromatinInGenes:
                        openChromatinInGenes[gene] += [loci]
                    else:
                        openChromatinInGenes[gene] = [loci]
                # Right border of window > Right border of open chromatin
                elif (tss[gene][1] + shift > tupel[1]) and (tss[gene][1] - shift < tupel[0]):
                    # Left border of window <= Left border of open chromatin ==> Window is larger than open chromatin
                    if gene in openChromatinInGenes:
                        openChromatinInGenes[gene] += [loci]
                    else:
                        openChromatinInGenes[gene] = [loci]
                # Left border of window <= Right border of open chromatin ==> Window enters open chromatin in the 5'
                # end stops in the tss window
                elif (tss[gene][1] + shift > tupel[1]) and (tss[gene][1] - shift >= tupel[0]) and (
                                tss[gene][1] - shift < tupel[1]):
                    if gene in openChromatinInGenes:
                        openChromatinInGenes[gene] += [loci]
                    else:
                        openChromatinInGenes[gene] = [loci]

    # Extract bound transcription factors
    affinities = extractTF_Affinity(openChromatinInGenes, args.affinity[0], tss, decay, loopsactivated, loopOCregions,
                                    geneloops, looplookuptable, loopdecay, loopcountscaling, countersiteonly,
                                    doublefeatures, enhancer)
    if decay:
        createAffinityFile(affinities, tfNames,
                           args.geneViewAffinity.replace("_Affinity_Gene_View.txt", "_Decay_Affinity_Gene_View.txt"),
                           tss, loopsactivated, doublefeatures)
        if args.sparseRep != 0:
            # noinspection PyPep8
            createSparseFile(affinities[0], tfNames, args.geneViewAffinity.replace("_Affinity_Gene_View.txt",
                                                                                   "_Decay_Sparse_Affinity_Gene_View.txt"),
                             tss, args.sparseRep)
    else:
        createAffinityFile(affinities, tfNames, args.geneViewAffinity, tss, loopsactivated, doublefeatures)
        if args.sparseRep != 0:
            createSparseFile(affinities[0], tfNames,
                             args.geneViewAffinity.replace("_Affinity_Gene_View.txt", "_Sparse_Affinity_Gene_View.txt"),
                             tss, args.sparseRep)

    if args.signalScale != "":
        scaledAffinities = extractTF_Affinity(openChromatinInGenes, args.signalScale, tss, decay, loopsactivated,
                                              loopOCregions, geneloops, looplookuptable, loopdecay, loopcountscaling,
                                              countersiteonly, doublefeatures, enhancer)
        if decay:
            # noinspection PyPep8
            createAffinityFile(scaledAffinities, tfNames, args.geneViewAffinity.replace("_Affinity_Gene_View.txt",
                                                                                        "_Decay_Scaled_Affinity_Gene_View.txt"),
                               tss, loopsactivated, doublefeatures)
            if args.sparseRep != 0:
                # noinspection PyPep8
                createSparseFile(scaledAffinities[0], tfNames, args.geneViewAffinity.replace("_Affinity_Gene_View.txt",
                                                                                             "_Decay_Scaled_Sparse_Affinity_Gene_View.txt"),
                                 tss, args.sparseRep)
        else:
            # noinspection PyPep8
            createAffinityFile(scaledAffinities, tfNames, args.geneViewAffinity.replace("_Affinity_Gene_View.txt",
                                                                                        "_Scaled_Affinity_Gene_View.txt"),
                               tss, loopsactivated, doublefeatures)
            if args.sparseRep != 0:
                # noinspection PyPep8
                createSparseFile(scaledAffinities[0], tfNames, args.geneViewAffinity.replace("_Affinity_Gene_View.txt",
                                                                                             "_Sparse_Scaled_Affinity_Gene_View.txt"),
                                 tss, args.sparseRep)

    if loopsactivated:
        createGenesWithLoopsFile(geneloops,
                                 args.geneViewAffinity.replace("_Affinity_Gene_View.txt", "_GenesWithLoops.txt"))

    now = datetime.datetime.now()
    print 'Finished annotation: ' + now.strftime("%Y-%m-%d-%H-%M-%S")


main()
