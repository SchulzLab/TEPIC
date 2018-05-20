import argparse

from utils import readIntraLoops, filterLoops, detectAllResolutions


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


def print_results(sample, loopOCregions, loops, resolution):
    loopcount = 0
    for chrKey in loops:
        loopcount += len(loops[chrKey])

    intersections = 0
    for intersection, sites in loopOCregions.items():
        if len(sites[0]) > 0 and len(sites[1]) > 0:
            intersections += 1

    print("{}\t{}\t{}\t{}\t{}".format(sample, resolution, intersections, round(float(intersections) / float(loopcount) * 100, 1), loopcount))


def main():
    parser = argparse.ArgumentParser(prog="calculateIntersectionRatio.py")
    parser.add_argument("regions", nargs=1, help="Open chromatin regions")
    parser.add_argument("loopfile", nargs=1,
                        help="If the name of the Hi-C loop file is provided, all open chromatin regions will be intersected with loop regions around the TSS of each gene.")
    parser.add_argument("--sample", type=str, default="-", help="Define a sample name to use in output.")
    args = parser.parse_args()

    oC = readOC_Region(args.regions[0])

    loops = readIntraLoops(args.loopfile[0])
    resolutions = detectAllResolutions(args.loopfile[0])

    print("Cell-line\thi-c_resolution\tintersection_count\tintersection_percentage\ttotal_count")

    for res in sorted(list(resolutions)):

        loops_of_res = dict(loops)
        filterLoops(loops_of_res, res)

        loopOCregions = intersectRegions(oC, loops_of_res)
        print_results(args.sample, loopOCregions, loops_of_res, res)

    loopOCregions = intersectRegions(oC, loops)
    # Repeat it once using loops from all resolutions
    print_results(args.sample, loopOCregions, loops, "All")


main()
