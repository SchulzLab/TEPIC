import argparse
import math
from operator import itemgetter

from SortedCollection import SortedCollection

# DEFINE IMPORTANT VARIABLES
SORTING_KEY = 1
WEIGHT_DECAY = 5000.0


def read_GTF(filename):
    gtf = open(filename, "r")
    tss = {}
    for l in gtf:
        s = l.split()
        if len(s) >= 9:
            if s[2] == "gene":
                chromosome = s[0].replace("chr", "")
                if chromosome not in tss:
                    tss[chromosome] = []
                if s[6] == "+":
                    tss[chromosome] += [(s[9], int(s[3]))]
                else:
                    tss[chromosome] += [(s[9], int(s[4]))]
    gtf.close()
    return tss


def read_GTF_Collection(filename):
    tss = {}
    with open(filename) as gtf_file:
        for l in gtf_file:
            s = l.split()
            if len(s) >= 9:
                if s[2] == "gene":
                    chromosome = s[0].replace("chr", "")
                    if chromosome not in tss:
                        tss[chromosome] = SortedCollection(key=itemgetter(SORTING_KEY))
                    if s[6] == "+":
                        tss[chromosome].insert_right((s[9], int(s[3])))
                    else:
                        tss[chromosome].insert_right((s[9], int(s[4])))
    return tss


def get_Gene_To_Chromosome_Mapping(filename):
    mapping = {}
    with open(filename) as gtf_file:
        for line in gtf_file:
            s = line.split()
            if len(s) >= 9:
                if s[2] == "gene":
                    mapping[s[9]] = s[0].replace("chr", "")
    return mapping


def read_Regions(filename):
    """
    :rtype: dictionary
    :param filename: string
    :return: regions as dictionary
    """
    regions = {}
    with open(filename) as region_file:
        rid = 0
        for entry in region_file:
            slr = entry.split()
            rchr = slr[0].replace("chr", "")
            rstart = int(slr[1])
            rend = int(slr[2])
            if rchr not in regions:
                regions[rchr] = []
            regions[rchr] += [(rid, rstart, rend)]
            rid += 1
        for chromosome in regions:
            regions[chromosome].sort()
    return regions


def read_Regions_Collection(filename):
    regions = {}
    with open(filename) as region_file:
        rid = 0
        for l in region_file:
            s = l.split()
            if len(s) >= 3:
                chromosome = s[0].replace("chr", "")
                rstart = int(s[1])
                rend = int(s[2])
                if chromosome not in regions:
                    regions[chromosome] = SortedCollection(key=itemgetter(SORTING_KEY))
                regions[chromosome].insert_right((rid, rstart, rend))
                rid += 1

    return regions


# Returns dictionary: {chrKey} -> [(start, end, coverage)]
def read_PeakCoverage(filename):
    peakcoverage = {}
    with open(filename) as peakcoverage_file:
        for line in peakcoverage_file:
            line = line.split()
            if len(line) < 4:
                raise ValueError("Error in peak coverage file! Expected: 'chr   start   end coverage' per line.")
            if line[0] not in peakcoverage:
                peakcoverage[line[0]] = []
            peakcoverage[line[0]] += [(int(line[1]), int(line[2]), float(line[3]))]
    return peakcoverage


def read_PeakCoverage_Collection(filename):
    peakcoverage = {}
    with open(filename) as peakcoverage_file:
        for line in peakcoverage_file:
            line = line.split()
            if len(line) < 4:
                raise ValueError("Error in peak coverage file! Expected: 'chr   start   end coverage' per line.")
            if line[0] not in peakcoverage:
                peakcoverage[line[0]] = SortedCollection(key=itemgetter(0, 1))
            peakcoverage[line[0]].insert_right((int(line[1]), int(line[2]), float(line[3])))
    return peakcoverage


def sanityCheck(a, b):
    # Perform some sanity checks
    for key in b:
        if key not in a:
            a[key] = SortedCollection()


def getNextMinima(collection, curr_minima, position, le):
    try:
        if le:
            next_gene = collection.find_le_index(position)
        else:
            next_gene = collection.find_gt_index(position)
    except ValueError:
        pass
    else:
        if curr_minima is not None:
            if abs(collection[curr_minima][SORTING_KEY] - position) > abs(collection[next_gene][SORTING_KEY] - position):
                return next_gene
        else:
            return next_gene
    return curr_minima


def getNearestNeighbors(annotations_collection, regions):
    gene_regions = {}
    for chromosome, chr_regions in regions.iteritems():
        collection = annotations_collection[chromosome]
        for region in chr_regions:
            minima = None
            position = region[1]
            minima = getNextMinima(collection, minima, position, True)
            minima = getNextMinima(collection, minima, position, False)
            position = region[2]
            minima = getNextMinima(collection, minima, position, True)
            minima = getNextMinima(collection, minima, position, False)

            if minima is not None:
                closest_gene = collection[minima]
                if closest_gene not in gene_regions:
                    gene_regions[closest_gene] = SortedCollection(key=itemgetter(SORTING_KEY))
                gene_regions[closest_gene].insert_right(region)

    return gene_regions


def getRegionsInWindow(annotations, regions_collection, window_size):
    gene_regions = {}
    for chromosome, chr_annotations in annotations.iteritems():
        collection = regions_collection[chromosome]
        for annotation in chr_annotations:
            gene_regions[annotation] = SortedCollection(key=itemgetter(SORTING_KEY))
            # Compute target interval indices
            lower_position = annotation[1] - window_size
            upper_position = annotation[1] + window_size
            try:
                left_item = collection.find_ge(lower_position)
            except ValueError:
                left_item = None
            try:
                right_item = collection.find_le(upper_position)
            except ValueError:
                right_item = None
            # Check if target interval is valid
            if left_item is not None and right_item is not None:
                left_index = collection.index(left_item)
                right_index = collection.index(right_item)
                if left_index <= right_index:
                    # Copy regions in target interval
                    for i in range(left_index, right_index + 1):
                        gene_regions[annotation].insert_right(collection[i])

    return gene_regions


def filterGeneRegions(gene_regions, gene_functions):
    filtered_gene_regions = {}
    for annotation, functional_regions in gene_functions.iteritems():
        regions = gene_regions[annotation]
        filtered_gene_regions[annotation] = SortedCollection(key=itemgetter(SORTING_KEY))
        for functional_region in functional_regions:
            try:
                left_boundary = regions.find_lt_index(functional_region[1])
            except ValueError:
                try:
                    left_boundary = regions.find_ge_index(functional_region[1])
                except ValueError:
                    left_boundary = len(regions)
            else:
                if regions[left_boundary][2] < functional_region[1]:
                    left_boundary += 1

            curr_index = left_boundary
            while curr_index < len(regions) and regions[curr_index][1] < functional_region[2]:
                filtered_gene_regions[annotation].insert_right(regions[curr_index])
                curr_index += 1
    return filtered_gene_regions


def computeDNaseFeatures(gene_regions, peakcoverage_collection, gene_to_chromosome, decay):
    gene_feature_matrix = {}
    weight_factor = 1.0
    for annotation, regions in gene_regions.iteritems():
        chromosome = gene_to_chromosome[annotation[0]]
        collection = peakcoverage_collection[chromosome]
        region_count = len(regions)
        total_length = 0.0
        total_signal = 0.0
        for region in regions:
            if decay:
                region_mid = int(((float(region[2]) - float(region[1])) / 2) + float(region[1]))
                weight_factor = math.exp(-(float(float(abs(annotation[1] - region_mid)) / WEIGHT_DECAY)))
            total_length += float(abs(region[2] - region[1])) * weight_factor
            try:
                cov_region = collection.find((region[1], region[2]))
            except ValueError:
                print "Skipping a region not found in peak coverage file..."
            else:
                total_signal += cov_region[2] * weight_factor
        gene_feature_matrix[annotation] = [region_count, total_length, total_signal]

    return gene_feature_matrix


def writeGeneFeatureMatrix(filename, feature_matrix):
    with open(filename, "w") as outputfile:
        outputfile.write("geneID    region_count    total_length    total_signal")
        outputfile.write('\n')
        for annotation, features in feature_matrix.iteritems():
            gene_id = annotation[0].split(".")[0].replace('"', '')
            outputfile.write(gene_id)
            for feature in features:
                outputfile.write('\t')
                outputfile.write(str(feature))
            outputfile.write('\n')


def main():
    parser = argparse.ArgumentParser(prog="computeRegionFeatures.py")
    parser.add_argument("gtf", help="Genome annotation file")
    parser.add_argument("regions", help="Region file")
    parser.add_argument("peakcoverage", help="Region coverage file.")
    parser.add_argument("--functional_regions", nargs="?", help="Enhancer/Promoter Region file")
    parser.add_argument("--window_size", nargs="?", type=int, help="Indicates whether to use a window of the given size for gene-region assignment. If no size is given, a nearest neighbor assignment is computed.")
    parser.add_argument("--outputprefix", nargs="?", default="", const="", help="Prefix to be used for output file(s).")
    decay_parser = parser.add_mutually_exclusive_group(required=True)
    decay_parser.add_argument('--decay', dest='decay', action='store_true', help="Flag option. If set an exponential decay function is used to weight the features. Default behavior is False.")
    decay_parser.add_argument('--no-decay', dest='decay', action='store_false')
    parser.set_defaults(decay=True)
    args = parser.parse_args()

    print "Loading files for DNaseFeatures..."
    # Extract TSS of GTF file
    annotations = read_GTF(args.gtf)
    annotations_collection = read_GTF_Collection(args.gtf)
    gene_to_chromosome = get_Gene_To_Chromosome_Mapping(args.gtf)

    # Load regions file
    regions = read_Regions(args.regions)
    regions_collection = read_Regions_Collection(args.regions)

    # Load Peak Coverage file
    peakcoverage_collection = read_PeakCoverage_Collection(args.peakcoverage)

    sanityCheck(annotations_collection, regions)
    sanityCheck(regions_collection, annotations)
    sanityCheck(peakcoverage_collection, annotations)

    print str(regions_collection.keys())
    print str(peakcoverage_collection.keys())

    print "Ready!"
    print "Starting annotation..."
    print "Using window size: " + str(args.window_size)
    print "Decay set to: " + str(args.decay)

    if args.window_size is not None:
        gene_regions = getRegionsInWindow(annotations, regions_collection, args.window_size)
        print "Assigned regions to " + str(len(gene_regions)) + " genes"
        count = 0
        for _, regs in gene_regions.iteritems():
            count += len(regs)
        print "Total number of assigned regions: " + str(count)
    else:
        gene_regions = getNearestNeighbors(annotations_collection, regions)
        print "Assigned regions to " + str(len(gene_regions)) + " genes"
        count = 0
        for _, regs in gene_regions.iteritems():
            count += len(regs)
        print "Total number of assigned regions: " + str(count)

    if args.functional_regions is not None:
        print "Loading functional regions..."
        if args.window_size is not None:
            functional_regions_collection = read_Regions_Collection(args.functional_regions)
            sanityCheck(functional_regions_collection, annotations)
            gene_functions = getRegionsInWindow(annotations, functional_regions_collection, args.window_size)
        else:
            functional_regions = read_Regions(args.functional_regions)
            sanityCheck(annotations_collection, functional_regions)
            gene_functions = getNearestNeighbors(annotations_collection, functional_regions)
        sanityCheck(gene_regions, gene_functions)
        print "Filtering regions using functional annotations..."
        gene_regions = filterGeneRegions(gene_regions, gene_functions)
        print "Assigned regions to " + str(len(gene_regions)) + " genes"
        count = 0
        for _, regs in gene_regions.iteritems():
            count += len(regs)
        print "Total number of assigned regions: " + str(count)

    feature_matrix = computeDNaseFeatures(gene_regions, peakcoverage_collection, gene_to_chromosome, args.decay)
    writeGeneFeatureMatrix(args.outputprefix + "RegionFeatures.txt", feature_matrix)

    print "Finished annotation!"

main()
