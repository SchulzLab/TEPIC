import argparse
import utils

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

def intersectRegions(oc, loops):
    looptable = {}
    for chrKey in oc:
        for octupel in oc[chrKey]:
            if (loops.has_key(chrKey)):
                chrLoops = loops[chrKey]
                for loop in chrLoops:
                    # open chromatin region in...
                    #	[  ----  ]
                    if (octupel[0] >= loop[1] and octupel[1] <= loop[2]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # --[------  ]
                    elif (octupel[0] <= loop[1] and octupel[1] <= loop[2] and octupel[1] >= loop[1]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # [  ------]--
                    elif (octupel[0] >= loop[1] and octupel[0] <= loop[2] and octupel[1] >= loop[2]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # --[--------]--
                    elif (octupel[0] < loop[1] and octupel[1] > loop[2]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][0].append(octupel)
                    # same procedure for second, counterpart loop-site
                    #	[  ----  ]
                    elif (octupel[0] >= loop[3] and octupel[1] <= loop[4]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # --[------  ]
                    elif (octupel[0] <= loop[3] and octupel[1] <= loop[4] and octupel[1] >= loop[3]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # [  ------]--
                    elif (octupel[0] >= loop[3] and octupel[0] <= loop[4] and octupel[1] >= loop[4]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)
                    # --[--------]--
                    elif (octupel[0] < loop[3] and octupel[1] > loop[4]):
                        if (not looptable.has_key(loop[0])):
                            looptable[loop[0]] = ([], [])
                        looptable[loop[0]][1].append(octupel)

    return looptable

def main():
    parser=argparse.ArgumentParser(prog="calculateIntersectionRatio.py")
    parser.add_argument("affinity",nargs=1,help="TRAP generated TF Affinity file")
    parser.add_argument("loopfile",nargs=1,help="If the name of the Hi-C loop file is provided, all open chromatin regions will be intersected with loop regions around the TSS of each gene.")
    args=parser.parse_args()

    oC=readOC_Region(args.affinity[0])

    loops = utils.readIntraLoops(args.loopfile[0])

    loopOCregions = intersectRegions(oC, loops)

    loopcount = 0
    for chrKey in loops:
        loopcount += len(loops[chrKey])

    intersections = 0
    for intersection in loopOCregions:
        sites = loopOCregions[intersection]
        if len(sites[0]) > 0 and len(sites[1]) > 0:
            intersections += 1

    print "Parsed " + str(loopcount) + " loops."
    print "Found " + str(intersections) + " loops intersecting in both loop-sites with affinity regions."
    print "That are " + str('%.1f' % round(float(intersections) / float(loopcount) * 100, 1)) + '%'


main()