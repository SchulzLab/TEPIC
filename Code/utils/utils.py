

# Reads a gtf file and generates a dictionary (key:gene, item:(#chromosom,TSS))
def readGTF(filename):
    gtf = open(filename, "r")
    tss = {}
    for l in gtf:
        s = l.split()
        if len(s) >= 9:
            if s[2] == "gene":
                if s[6] == "+":
                    tss[s[9]] = (s[0].replace("chr", ""), int(s[3]))
                else:
                    tss[s[9]] = (s[0].replace("chr", ""), int(s[4]))
    gtf.close()
    return tss


# Reads a loop file and returns a dictionary containing all intrachromosomal loops per chromosome in a list:
# {key : value} -> {chr : [(loopID, start X, end X, start Y, end Y, observations count, resolution)]}
def readIntraLoops(loopsFile):
    lf = open(loopsFile, 'r')
    loops = {}

    loopID = 0
    for l in lf:
        s = l.split()
        if len(s) >= 8:
            if s[0] == s[3]:  # check if loop is intra-chromosomal
                loopID += 1
                resolution = int(s[2]) - int(s[1])
                if s[0] not in loops:
                    loops[s[0]] = []
                loops[s[0]].append((loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7]), resolution))
    lf.close()
    return loops


# Reads a loop file and returns a dictionary containing all intrachromosomal loops accessible by an unique ID:
# {key : value} -> {loopID : (loopID, start X, end X, start Y, end Y, observations count, resolution, chromosome)}
# we keep the loop format for consistent programming and accessing tupel fields
def readIntraLoopsLookupTable(loopsFile):
    lf = open(loopsFile, 'r')
    loops = {}

    loopID = 0
    for l in lf:
        s = l.split()
        if len(s) >= 8:
            if s[0] == s[3]:  # check if loop is intra-chromosomal
                loopID += 1
                resolution = int(s[2]) - int(s[1])
                loops[loopID] = (loopID, int(s[1]), int(s[2]), int(s[4]), int(s[5]), int(s[7]), resolution, s[0])
    lf.close()
    return loops


def readEnhancer(enhancerFile):
    ef = open(enhancerFile, 'r')
    enhancer = {}

    for l in ef:
        s = l.split()
        if len(s) >= 3:
            s[0] = s[0].replace("chr", "")
            if s[0] not in enhancer:
                enhancer[s[0]] = []
            enhancer[s[0]].append((int(s[1]), int(s[2])))
    ef.close()
    return enhancer


def detectAllResolutions(loopsFile):
    lf = open(loopsFile, 'r')
    resolutions = set()
    lf.readline()

    for l in lf:
        s = l.split()
        if len(s) >= 8:
            resolution = int(s[2]) - int(s[1])
            resolutions.add(resolution)
    lf.close()
    return resolutions


# Reads a loop file and returns a dictionary containing all intrachromosomal loops per chromosome in a list:
# TODO: define format
def readInterLoops(loopsFile):
    # TODO: implement this
    return 0


# Writes given header and body to a file defined by filename
def writeToFile(filename, header, body):
    f = open(filename, 'w')
    f.write(str(header))
    f.write('\n')
    f.write(str(body))
    f.close()


def filterLoops(loops, resolution):
    for chrKey in loops:
        chrLoops = loops[chrKey]
        for loop in chrLoops:
            if (loop[2] - loop[1]) != resolution:
                chrLoops.remove(loop)
