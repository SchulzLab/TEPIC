import argparse
import math
from operator import itemgetter, add

from SortedCollection import SortedCollection

# DEFINE IMPORTANT VARIABLES
SORTING_KEY = 1
SORTING_KEY_PAIRED = 3
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
            chromosome = line[0].replace("chr", "")
            if chromosome not in peakcoverage:
                peakcoverage[chromosome] = []
            peakcoverage[chromosome] += [(int(line[1]), int(line[2]), float(line[3]))]
    return peakcoverage


def read_PeakCoverage_Collection(filename):
    peakcoverage = {}
    with open(filename) as peakcoverage_file:
        for line in peakcoverage_file:
            line = line.split()
            if len(line) < 4:
                raise ValueError("Error in peak coverage file! Expected: 'chr   start   end coverage' per line.")
            chromosome = line[0].replace("chr", "")
            if chromosome not in peakcoverage:
                peakcoverage[chromosome] = SortedCollection(key=itemgetter(0, 1))
            peakcoverage[chromosome].insert_right((int(line[1]), int(line[2]), float(line[3])))
    return peakcoverage


def read_Intra_Loop_Regions_Collection(filename, resolution_filter):
    left_region_collection = {}
    right_region_collection = {}
    with open(filename) as hi_c_file:
        # skip header
        hi_c_file.readline()
        loopID = 0
        for line in hi_c_file:
            line = line.split()
            if len(line) >= 8:
                if line[0] == line[3]:  # check if loop is intra-chromosomal
                    loopID += 1
                    chromosome = line[0]
                    resolution = int(line[2]) - int(line[1])
                    if chromosome not in left_region_collection:
                        left_region_collection[chromosome] = SortedCollection(key=itemgetter(SORTING_KEY))
                        right_region_collection[chromosome] = SortedCollection(key=itemgetter(SORTING_KEY))
                    if resolution_filter is None or resolution_filter == resolution:
                        left_region_collection[chromosome].insert_right((loopID, int(line[1]), int(line[2]), int(line[4]), int(line[5])))
                        right_region_collection[chromosome].insert_right((loopID, int(line[4]), int(line[5]), int(line[1]), int(line[2])))
    return left_region_collection, right_region_collection


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
            left_minima = getNextMinima(collection, None, region[1], True)
            left_minima = getNextMinima(collection, left_minima, region[1], False)

            right_minima = getNextMinima(collection, None, region[2], True)
            right_minima = getNextMinima(collection, right_minima, region[2], False)

            closest_gene = None
            if left_minima is not None and right_minima is not None:
                if abs(collection[left_minima][SORTING_KEY] - region[1]) > abs(collection[right_minima][SORTING_KEY] - region[2]):
                    closest_gene = collection[right_minima]
                else:
                    closest_gene = collection[left_minima]
            elif left_minima is not None:
                closest_gene = collection[left_minima]
            elif right_minima is not None:
                closest_gene = collection[right_minima]

            if closest_gene is not None:
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
                left_item = collection.find_lt(lower_position)
            except ValueError:
                try:
                    left_item = collection.find_ge(lower_position)
                except ValueError:
                    left_item = None
            else:
                if left_item[SORTING_KEY + 1] < lower_position:
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


def filterGeneRegions(regions_collection, gene_functions, gene_to_chromosome=None, paired_functions=False):
    filtered_gene_regions = {}
    filtered_gene_functions = {}
    for annotation, functional_regions in gene_functions.iteritems():
        if gene_to_chromosome is None:
            regions = regions_collection[annotation]
        else:
            regions = regions_collection[gene_to_chromosome[annotation[0]]]
        (filtered_gene_functions[annotation], filtered_gene_regions[annotation]) = get_intersecting_regions(functional_regions, regions, paired_a=paired_functions)
    return filtered_gene_regions, filtered_gene_functions


def get_intersecting_regions(a_regions, b_collection, paired_a=False):
    intersection_a = SortedCollection(key=itemgetter(SORTING_KEY))
    intersection_b = SortedCollection(key=itemgetter(SORTING_KEY))
    for a_region in a_regions:
        try:
            left_boundary = b_collection.find_lt_index(a_region[SORTING_KEY])
        except ValueError:
            try:
                left_boundary = b_collection.find_ge_index(a_region[SORTING_KEY])
            except ValueError:
                left_boundary = len(b_collection)
        else:
            if b_collection[left_boundary][SORTING_KEY+1] < a_region[SORTING_KEY]:
                left_boundary += 1

        curr_index = left_boundary
        if curr_index < len(b_collection) and b_collection[curr_index][SORTING_KEY] <= a_region[SORTING_KEY+1]:
            intersection_a.insert_right(a_region)
            while curr_index < len(b_collection) and b_collection[curr_index][SORTING_KEY] <= a_region[SORTING_KEY+1]:
                try:
                    intersection_b.find(b_collection[curr_index][SORTING_KEY])
                except ValueError:
                    intersection_b.insert_right(b_collection[curr_index])
                curr_index += 1

    if paired_a:
        for a_region in a_regions:
            try:
                left_boundary = b_collection.find_lt_index(a_region[SORTING_KEY_PAIRED])
            except ValueError:
                try:
                    left_boundary = b_collection.find_ge_index(a_region[SORTING_KEY_PAIRED])
                except ValueError:
                    left_boundary = len(b_collection)
            else:
                if b_collection[left_boundary][SORTING_KEY+1] < a_region[SORTING_KEY_PAIRED]:
                    left_boundary += 1

            curr_index = left_boundary
            if curr_index < len(b_collection) and b_collection[curr_index][SORTING_KEY] <= a_region[SORTING_KEY_PAIRED+1]:
                intersection_a.insert_right(a_region)
                while curr_index < len(b_collection) and b_collection[curr_index][SORTING_KEY] <= a_region[SORTING_KEY_PAIRED+1]:
                    try:
                        intersection_b.find(b_collection[curr_index][SORTING_KEY])
                    except ValueError:
                        intersection_b.insert_right(b_collection[curr_index])
                    curr_index += 1

    return intersection_a, intersection_b


def merge_gene_regions(gene_regions, add_gene_regions):
    for annotation, regions in add_gene_regions.iteritems():
        if annotation not in gene_regions:
            gene_regions[annotation] = SortedCollection(key=itemgetter(SORTING_KEY))
        for region in regions:
            try:
                gene_regions[annotation].find(region[SORTING_KEY])
            except ValueError:
                gene_regions[annotation].insert_right(region)
    return gene_regions


def remove_gene_regions(gene_regions, subs_gene_regions):
    for annotation, regions in subs_gene_regions.iteritems():
        if annotation in gene_regions:
            for region in regions:
                try:
                    gene_regions[annotation].remove(region)
                except ValueError:
                    pass
    return gene_regions


def computeRegionFeatures(gene_regions, peakcoverage_collection, gene_to_chromosome, decay):
    gene_feature_matrix = {}
    weight_factor = 1.0
    for annotation, regions in gene_regions.iteritems():
        chromosome = gene_to_chromosome[annotation[0]]
        collection = peakcoverage_collection[chromosome]
        region_count = 0.0
        total_length = 0.0
        total_signal = 0.0
        for region in regions:
            if decay:
                region_mid = int(((float(region[2]) - float(region[1])) / 2) + float(region[1]))
                weight_factor = math.exp(-(float(float(abs(annotation[1] - region_mid)) / WEIGHT_DECAY)))
            region_count += weight_factor
            total_length += float(abs(region[2] - region[1])) * weight_factor
            try:
                cov_region = collection.find((region[1], region[2]))
            except ValueError:
                print "Skipping a region not found in peak coverage file..."
            else:
                total_signal += cov_region[2] * weight_factor
        gene_feature_matrix[annotation] = [region_count, total_length, total_signal]

    return gene_feature_matrix


def mergeGeneFeatureMatrices(matrix_a, matrix_b):
    for annotation in matrix_b:
        if annotation in matrix_a:
            matrix_a[annotation] = map(add, matrix_a[annotation], matrix_b[annotation])
        else:
            matrix_a[annotation] = matrix_b[annotation]


def writeGeneFeatureMatrix(filename, feature_matrix, second_matrix=None):
    with open(filename, "w") as outputfile:
        outputfile.write("geneID\tregion_count\ttotal_length\ttotal_signal")
        if second_matrix is not None:
            outputfile.write("\tregion_count_loop\ttotal_length_loop\ttotal_signal_loop")
        outputfile.write('\n')
        for annotation, features in feature_matrix.iteritems():
            gene_id = annotation[0].split(".")[0].replace('"', '')
            outputfile.write(gene_id)
            for feature in features:
                outputfile.write('\t')
                outputfile.write(str(feature))
            if second_matrix is not None:
                for feature in second_matrix[annotation]:
                    outputfile.write('\t')
                    outputfile.write(str(feature))
            outputfile.write('\n')


def main():
    parser = argparse.ArgumentParser(prog="computeRegionFeatures.py")
    parser.add_argument("gtf", help="Genome annotation file")
    parser.add_argument("regions", help="Region file")
    parser.add_argument("peakcoverage", help="Region coverage file.")
    parser.add_argument("--functional_regions", nargs="?", default=None, help="Enhancer/Promoter Region file")
    parser.add_argument("--window_size", nargs="?", type=int, help="Indicates whether to use a window of the given size for gene-region assignment. If no size is given, a nearest neighbor assignment is computed.")
    parser.add_argument("--outputprefix", nargs="?", default="", const="", help="Prefix to be used for output file(s).")
    decay_parser = parser.add_mutually_exclusive_group(required=True)
    decay_parser.add_argument('--decay', dest='decay', action='store_true', help="Flag option. If set an exponential decay function is used to weight the features. Default behavior is False.")
    decay_parser.add_argument('--no-decay', dest='decay', action='store_false')
    parser.set_defaults(decay=False)
    parser.add_argument("--hi_c_regions", nargs="?", default=None, help="If the name of the Hi-C loop file is provided, all regions will be intersected with annotated Hi-C (intrachromosomal) loop regions around the TSS of each annotated gene.")
    parser.add_argument("--loopwindows", nargs="?", default=25000, type=int, help="Defines the window-size around the TSS in which all loops are considered for intersecting with openChromatin regions.")
    parser.add_argument("--resolution", nargs="?", default=None, help="Defines the Hi-C resolution of the loops which should be considered.")
    hic_decay_parser = parser.add_mutually_exclusive_group(required=False)
    hic_decay_parser.add_argument('--hi_c-decay', dest='hic_decay', action='store_true', help="Flag option. If set an exponential decay function is used to weight the Hi-C features. Default behavior is False.")
    hic_decay_parser.add_argument('--no_hi_c-decay', dest='hic_decay', action='store_false')
    parser.set_defaults(hic_decay=False)
    double_feature_parser = parser.add_mutually_exclusive_group(required=False)
    double_feature_parser.add_argument('--double_features', dest='double_features', action='store_true', help="Flag option. If set separate gene-wise features are computed for promoter and loop regions instead of aggregating all scores.")
    double_feature_parser.add_argument('--no_double_features', dest='double_features', action='store_false')
    parser.set_defaults(double_features=False)
    args = parser.parse_args()

    print "Loading files ..."
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
    print "Hi-C Decay set to: " + str(args.hic_decay)

    if args.window_size is not None:
        gene_regions = getRegionsInWindow(annotations, regions_collection, args.window_size/2)
    else:
        gene_regions = getNearestNeighbors(annotations_collection, regions)

    print "Assigned regions to " + str(len(gene_regions)) + " genes"
    count = 0
    for _, regs in gene_regions.iteritems():
        count += len(regs)
    print "Total number of assigned regions: " + str(count)

    gene_hic_feature_matrix = None
    if args.functional_regions is not None and args.hi_c_regions is None:
        print "Loading functional regions..."
        if args.window_size is not None:
            functional_regions_collection = read_Regions_Collection(args.functional_regions)
            sanityCheck(functional_regions_collection, annotations)
            gene_functions = getRegionsInWindow(annotations, functional_regions_collection, args.window_size/2)
        else:
            functional_regions = read_Regions(args.functional_regions)
            sanityCheck(annotations_collection, functional_regions)
            gene_functions = getNearestNeighbors(annotations_collection, functional_regions)
        sanityCheck(gene_regions, gene_functions)
        print "Filtering regions using functional annotations..."
        (gene_regions, _) = filterGeneRegions(gene_regions, gene_functions)

        print "Assigned regions to " + str(len(gene_regions)) + " genes"
        count = 0
        for _, regs in gene_regions.iteritems():
            count += len(regs)
        print "Total number of assigned regions: " + str(count)

    elif args.hi_c_regions is not None:
        print "Running into Hi-c mode..."
        print "Using loop-window size: " + str(args.loopwindows)
        # read hi-c file and keep loops with desired resolution
        resolution = None
        if args.resolution is not None:
            resolution = int(args.resolution)
        print "Using loop resolution " + str(resolution)
        (regions_left_collection, regions_right_collection) = read_Intra_Loop_Regions_Collection(args.hi_c_regions, resolution)
        sanityCheck(regions_left_collection, annotations)
        sanityCheck(regions_right_collection, annotations)
        # get gene-wise hi-c loops in loopwindows
        gene_loops_left = getRegionsInWindow(annotations, regions_left_collection, args.loopwindows/2)
        gene_loops_right = getRegionsInWindow(annotations, regions_right_collection, args.loopwindows/2)
        # filter regions that are in these gene-wise hi-c regions
        (gene_regions_left, _) = filterGeneRegions(regions_collection, gene_loops_left, gene_to_chromosome, paired_functions=True)
        (gene_regions_right, _) = filterGeneRegions(regions_collection, gene_loops_right, gene_to_chromosome, paired_functions=True)
        if args.functional_regions is not None:
            print "Using also enhancer regions..."
            # !use regions in tss window + intersect functional regions with hi-c regions and intersect those again with regions!
            # read functional regions
            functional_regions_collection = read_Regions_Collection(args.functional_regions)
            sanityCheck(functional_regions_collection, annotations)
            # filter the regions in gene-wise hi-c regions further using the functional regions
            (_, gene_regions_left) = filterGeneRegions(functional_regions_collection, gene_regions_left, gene_to_chromosome)
            (_, gene_regions_right) = filterGeneRegions(functional_regions_collection, gene_regions_right, gene_to_chromosome)
        # OR: !just use regions in tss window and those in hi-c regions!
        # merge gene regions from left and right Hi-C regions
        gene_hic_regions = merge_gene_regions(gene_regions_left, gene_regions_right)
        # Filter duplicated regions resulting from multiple region searches
        gene_hic_regions = remove_gene_regions(gene_hic_regions, gene_regions)

        print "Assigned regions to " + str(len(gene_hic_regions)) + " genes"
        count = 0
        for _, regs in gene_hic_regions.iteritems():
            count += len(regs)
        print "Total number of assigned regions: " + str(count)
        gene_hic_feature_matrix = computeRegionFeatures(gene_hic_regions, peakcoverage_collection, gene_to_chromosome, args.hic_decay)

    gene_feature_matrix = computeRegionFeatures(gene_regions, peakcoverage_collection, gene_to_chromosome, args.decay)
    if gene_hic_feature_matrix is not None and not args.double_features:
        mergeGeneFeatureMatrices(gene_feature_matrix, gene_hic_feature_matrix)
        gene_hic_feature_matrix = None

    output_filename = "RegionFeatures"
    if args.outputprefix:
        output_split = args.outputprefix.rsplit("/", 1)
        output_filename = args.outputprefix
        if len(output_split) > 1 and not output_split[1]:
            output_filename = args.outputprefix + "RegionFeatures"

    writeGeneFeatureMatrix(output_filename + ".txt", gene_feature_matrix, second_matrix=gene_hic_feature_matrix)

    print "Finished annotation!"

main()
