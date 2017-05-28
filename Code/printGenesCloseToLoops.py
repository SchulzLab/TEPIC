import argparse
import utils


def readGTF(filename):
    gtf = open(filename, "r")
    tss = {}
    for l in gtf:
        s = l.split()
        if len(s) >= 9:
            if s[2] == "gene":
                if s[6] == "+":
                    tss[s[9]] = (s[0].replace("chr", ""), int(s[3]), l.replace('\n', ''))
                else:
                    tss[s[9]] = (s[0].replace("chr", ""), int(s[4]), l.replace('\n', ''))
    gtf.close()
    return tss


def findLoopsNearbyGenes(tss, loops, window_size):
    loopwindows = int(window_size / 2)
    for geneID in tss:
        startpos = tss[geneID][1]
        if tss[geneID][0] in loops:
            chrLoops = loops[tss[geneID][0]]
            for loop in chrLoops:

                leftmiddle = int(((float(loop[2]) - float(loop[1])) / 2) + float(loop[1]))
                if abs(leftmiddle - startpos) <= loopwindows:
                    print tss[geneID][2]
                    break

                else:
                    rightmiddle = int(((float(loop[3]) - float(loop[4])) / 2) + float(loop[3]))
                    if abs(rightmiddle - startpos) <= loopwindows:
                        print tss[geneID][2]
                        break


def main():
    parser = argparse.ArgumentParser(prog="printGenesCloseToLoops.py")
    parser.add_argument("gtf", help="Genome annotation file")
    parser.add_argument("loopfile", help="If the name of the Hi-C loop file is provided, all open chromatin regions will be intersected with loop regions around the TSS of each gene. Activates the Hi-C mode.")
    parser.add_argument("window_size", help="Defines the window-size around the TSS in which all loops are considered for intersecting with openChromatin regions.", type=int)

    args = parser.parse_args()
    gtf = readGTF(args.gtf)
    loops = utils.readIntraLoops(args.loopfile)
    findLoopsNearbyGenes(gtf, loops, args.window_size)


main()
