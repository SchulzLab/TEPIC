import argparse
import copy

from utils import detectAllResolutions, filterLoops, readGTF, readIntraLoops, writeToFile


def run(genes, intraLoops, windows):
    matchedLoops = {}

    for geneKey in genes:
        chr = genes[geneKey][0].upper()
        if chr != "X" and chr != "Y" and chr != "M":
            chrLoops = intraLoops[chr]
            pos = genes[geneKey][1]
            if geneKey not in matchedLoops:
                matchedLoops[geneKey] = {}

            for loop in chrLoops:
                resolution = loop[6]

                for window in windows:
                    posRight = pos + window
                    posLeft = pos - window

                    if window not in matchedLoops[geneKey]:
                        matchedLoops[geneKey][window] = {}

                    if resolution not in matchedLoops[geneKey][window]:
                        matchedLoops[geneKey][window][resolution] = []

                    if posLeft <= loop[1] <= posRight or posLeft <= loop[2] <= posRight or (loop[1] <= posLeft and posRight <= loop[2]):
                        matchedLoops[geneKey][window][resolution].append(loop)

                    else:
                        if posLeft <= loop[3] <= posRight or posLeft <= loop[4] <= posRight or (loop[3] <= posLeft and posRight <= loop[4]):
                            matchedLoops[geneKey][window][resolution].append(loop)

    return matchedLoops


def postProcessing(results, sample, resolution):
    per_window_counter = {}

    for geneKey, windows in results.items():

        for windowKey, resolutions in windows.items():

            for resolutionKey, hits in resolutions.items():
                if resolution and resolutionKey != resolution:
                    raise ValueError("Resolution " + resolution + "cannot be found in the data!")
                else:
                    if len(hits):
                        if windowKey not in per_window_counter:
                            per_window_counter[windowKey] = 0
                        per_window_counter[windowKey] += 1

    if not resolution:
        resolution = "All"

    for window, gene_count in per_window_counter.items():
        output = str(sample) + "\t"
        output += str(resolution) + "\t"
        output += str(window) + '\t'
        output += str(gene_count)
        print(output)


def collectResolutionsPerWindowToSize(annotationFile, loopsFile, windows, sample):
    print('Indexing TSS')
    tss = readGTF(annotationFile)
    print('Found ' + str(len(tss)) + " genes")
    gene_count = 0
    for _, gene in tss.items():
        if gene[0].upper() != "X" and gene[0].upper() != "Y" and gene[0].upper() != "M":
            gene_count += 1
    print('Found ' + str(gene_count) + " genes on autosomes")
    print('Indexing Loops')
    loops = readIntraLoops(loopsFile)
    print('Preprocessing')
    for window in windows:
        if window < 100:
            print('Window radius too small... please use greater values e.g. 100 and above.')
            return 1

    resolutions = detectAllResolutions(loopsFile)
    print('Running core algorithm')
    print('sample\thi-c_resolution\twindow\tgene_count')
    for res in sorted(list(resolutions)):

        loops_of_res = copy.deepcopy(loops)
        filterLoops(loops_of_res, res)
        results = run(tss, loops_of_res, windows)
        postProcessing(results, sample, res)

    results = run(tss, loops, windows)
    postProcessing(results, sample, None)

    return 0


def convert_multipliers(n_str):
    last_char= n_str[-1:].upper()
    if last_char == 'K':
        return n_str[:-1]  + "000"
    elif last_char == 'M':
        return n_str[:-1] + "000000"
    else:
        return n_str


######
##
# #	Required arguments: 1) annotation File (.gtf), 2) loop file in Hi-C loop format,
# # 3) a comma separated enumeration of window radii (window radii around TSS in thousand bases)
##
######

# preparation-routine 
parser = argparse.ArgumentParser(
    description='Collects all loops in a window around genes and extracts some statistical values')
parser.add_argument('annotation', help='Path to an annotation file')
parser.add_argument('loops', help='Path to a loop file')
parser.add_argument('windows', type=str, help='A list of window radii around the genestart which should be scanned')
parser.add_argument("--sample", type=str, default="-", help="Define a sample name to use in output.")

args = parser.parse_args()
win = args.windows.split(',')
win = [int(convert_multipliers(numeric_string)) for numeric_string in win]

print('Starting to collect data...')
collectResolutionsPerWindowToSize(args.annotation, args.loops, win, args.sample)
print('\n-> Completed all!')
