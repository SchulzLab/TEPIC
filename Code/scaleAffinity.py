import sys
import argparse as argp
import traceback as trb


def parse_commandline():
    """
    """
    parser = argp.ArgumentParser(prog='scaleAffinity.py')
    parser.add_argument('-s', '--scales', type=str, dest='signalScales',
                        help='File containing scale information')
    parser.add_argument('-a', '--affinities', type=str, dest='affinity',
                        help='TAB-separated TF affinity file')
    parser.add_argument('--is-sorted', default=False, action='store_true', dest='issorted',
                        help='File with scale information is in same order as regions in affinty file')
    parser.add_argument('--scale-col', type=int, default=4, dest='scalecol',
                        help='1-based index of column containing scaling factor. Default: 4'
                             ' (i.e., simple chrom-start-end-value file is assumed)')
    args = parser.parse_args()
    return args


def scaling_sorted(args):
    """
    """
    scale_idx = args.scalecol - 1  # Python: 0-based indexing
    with open(args.affinity, 'r') as affn:
        header = affn.readline()
        sys.stdout.write(header)
        with open(args.signalScales, 'r') as scales:
            for aln, sln in zip(affn, scales):
                factor = float(sln.split()[scale_idx])
                annreg = aln.strip().split()
                scaled = [str(float(x) * factor) for x in annreg[1:]]
                sys.stdout.write(annreg[0] + '\t' + '\t'.join(scaled) + '\n')
    return


def check_header(fpath):
    """
    :param fpath:
    :return:
    """
    has_header = False
    with open(fpath, 'r') as bedlike:
        hd = bedlike.readline().strip()
        if hd.startswith('#'):
            has_header = True
        else:
            try:
                cols = hd.split()
                _ = int(cols[1])
                _ = int(cols[2])
            except ValueError:
                has_header = True
    return has_header


def scaling_lookup(args):
    """
    """
    assert args.scalecol > 3, 'Unsorted input requires the scales file to be BED-like, i.e., ' \
                              'chrom(1) - start(2) - end(3) have to be the first three columns; ' \
                              'yet you provided {} as column index for the scaling factor.'.format(args.scalecol)
    skip_header = check_header(args.signalScales)
    coord_lut = dict()
    scale_idx = args.scalecol - 1  # Python: 0-based indexing
    with open(args.signalScales, 'r') as scales:
        if skip_header:
            _ = scales.readline()
        for line in scales:
            try:
                cols = line.split()
                k = cols[0] + ':' + cols[1] + '-' + cols[2]
            except IndexError:
                # skip over empty lines
                continue
            # notably, if non-numeric column was specified,
            # this raises
            coord_lut[k] = float(cols[scale_idx])
    with open(args.affinity, 'r') as affn:
        header = affn.readline()
        sys.stdout.write(header)
        for line in affn:
            annreg = line.strip().split()
            factor = float(coord_lut[annreg[0]])
            scaled = [str(float(x) * factor) for x in annreg[1:]]
            sys.stdout.write(annreg[0] + '\t' + '\t'.join(scaled) + '\n')
    return


def main():
    """
    """
    args = parse_commandline()
    assert args.scalecol >= 0, 'Index for column containing scaling factor has to be >= 0; ' \
                               'you specified {}'.format(args.scalecol)
    if args.issorted:
        scaling_sorted(args)
    else:
        scaling_lookup(args)
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\n{}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)

