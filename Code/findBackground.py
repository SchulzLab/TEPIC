# coding=utf-8

import sys as sys
import argparse as argp
import collections as col
import csv as csv
import functools as fnt
import io as io
import multiprocessing as mp
import operator as op
import threading as thd
import traceback as trb
import time as ti

# imports of non-standard libraries
# this represents the set of external
# dependencies for this script
import numpy as np
import numpy.random as rng
import scipy.spatial as spat
import twobitreader as tbr

"""
Stand-alone script that expects a BED-like input file
and searches for a similar genomic region for each region
in the input file. In other words, given a list of regions
of interest (the foreground), the scripts searches for a
similar set of background regions.
Similarity is defined - for this script - in terms of GC
content and region length.

Note that this code has been extracted from another project
and has been designed for a slightly different use case.
This implies that performance is not optimized for the
specific use case in TEPIC.
"""

__DEVELOPER__ = 'Peter Ebert'
__DEV_EMAIL__ = 'pebert@mpi-inf.mpg.de'
__MAINTAINER__ = 'Florian Schmidt'
__MAINT_EMAIL__ = 'fschmidt@mmci.uni-saarland.de'


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog='TEPIC/findBackground.py', add_help=True)

    ag = parser.add_argument_group('Input/output parameters')
    ag.add_argument('--input', '-i', type=str, dest='inputfile', required=True,
                    help='Path to input file. First three columns in file'
                         ' are expected to be chrom - start - end.')
    ag.add_argument('--genome', '-g', type=str, dest='genome', required=True,
                    help='Path to genome reference file in 2bit format.')
    ag.add_argument('--output', '-o', type=str, dest='outputfile', default='stdout',
                    help='Path to output file or stdout. Default: stdout')

    ag = parser.add_argument_group('Runtime parameters')
    ag.add_argument('--workers', '-w', type=int, dest='workers', default=1,
                    help='Number of CPU cores to use. 1 CPU core'
                         ' processes 1 chromosome at a time. Default: 1')
    ag.add_argument('--time-out', '-to', type=int, dest='timeout', default=3,
                    help='Maximal number of minutes to spend searching for'
                         ' background regions per chromosome. Default: 3 minutes')
    ag.add_argument('--threshold', '-th', type=int, dest='threshold', default=90,
                    help='Stop searching after having found more than <THRESHOLD>%%'
                         ' matches per chromosome. Default: 90%%')
    ag.add_argument('--eps-init', '-ei', type=float, dest='epsinit', default=1.,
                    help='Init value for epsilon. Error tolerance in percentage points'
                         ' for similarity matching. Default: 1.0 ppt')
    ag.add_argument('--eps-step', '-es', type=float, dest='epsstep', default=0.5,
                    help='Increment epsilon at each iteration by this value. Default: 0.5')
    ag.add_argument('--eps-max', '-em', type=float, dest='epsmax', default=2.,
                    help='Maximal value for epsilon. After reaching this value, restart'
                         ' search with different starting positions. Default: 2.0')
    return parser.parse_args()


def read_input_file(fpath):
    """
    :param fpath:
    :return:
    """
    foreground = col.defaultdict(list)
    with open(fpath, 'r') as inputfile:
        rows = csv.reader(inputfile, delimiter='\t')
        try:
            for idx, row in enumerate(rows, start=1):
                region = {'name': idx, 'chrom': row[0],
                          'start': int(row[1]), 'end': int(row[2]),
                          'length': int(row[2]) - int(row[1])}
                foreground[row[0]].append(region)
        except (ValueError, IndexError, TypeError) as err:
            sys.stderr.write('\nError processing line {} in input file {}\nRecord: {}\n'.format(idx, fpath, row))
            raise err
    assert foreground, 'No input regions read from file {}'.format(fpath)
    return foreground


def load_chromosome_sequence(chrom, fpath):
    """
    :param chrom:
    :param fpath:
    :return:
    """
    tbf = tbr.TwoBitFile(fpath)
    for key in tbf.keys():
        tbf[key.replace("chr","")]=tbf[key]
    assert chrom in tbf, 'Chromosome {} not in reference file {}'.format(chrom, fpath)
    chrom_seq = tbf[chrom]
    return str(chrom_seq).upper()


def compute_data_limits(regions, rt_params):
    """
    :param regions:
    :param rt_params:
    :return:
    """
    fac_hi = 1 + rt_params['epsmax'] / 100.
    fac_lo = 1 - rt_params['epsmax'] / 100.
    reglen_max = max([r['length'] for r in regions])
    reglen_min = min([r['length'] for r in regions])
    len_bound_hi = np.ceil(reglen_max * fac_hi)
    len_bound_lo = np.floor(reglen_min * fac_lo)
    threshold = min(len(regions), np.ceil(len(regions) * (rt_params['threshold'] / 100.)))
    limits = {'factor_lo': fac_lo, 'factor_hi': fac_hi,
              'len_bound_lo': len_bound_lo, 'len_bound_hi': len_bound_hi,
              'threshold': threshold, 'num_samples': len(regions)}
    return limits


def bulk_add_seq_features(regions, sequence, yardstick):
    """
    :param regions:
    :param sequence:
    :param yardstick:
    :return:
    """
    regions = sorted(regions, key=lambda x: x['name'])
    idxmap = dict()
    features = []
    taken = np.zeros(len(sequence), dtype=np.bool)
    for row, rec in enumerate(regions):
        ft_pct_len = np.round(rec['length'] / yardstick * 100., 2)
        seq = sequence[rec['start']:rec['end']]
        num_c = seq.count('C')
        num_g = seq.count('G')
        ft_pct_gc = np.round((num_g + num_c) / float(rec['length']) * 100., 2)
        idxmap[row] = rec['name']
        features.append([ft_pct_len, ft_pct_gc])
        taken[rec['start']:rec['end']] = 1
    features = np.array(features, dtype=np.float32)
    return features, taken, idxmap


def compute_seq_features(yardstick, sequence):
    """
    :param yardstick:
    :param sequence:
    :return:
    """
    seq_len = float(len(sequence))
    ft_pct_len = np.round(seq_len / yardstick * 100., 2)
    num_c = sequence.count('C')
    num_g = sequence.count('G')
    ft_pct_gc = np.round((num_g + num_c) / seq_len * 100., 2)
    return [ft_pct_len, ft_pct_gc]


def start_relaxed_search(nntree, chromseq, covered, idxmap, limits, params):
    """
    :param nntree:
    :param chromseq:
    :param covered:
    :param idxmap:
    :param limits:
    :param params:
    :return:
    """
    timeout = thd.Event()
    if params['timeout'] > 0:
        timer = thd.Timer(params['timeout'] * 60, timeout.set)
        timer.start()
    chrom_bound = len(chromseq)
    coordinates = np.arange(chrom_bound, dtype=np.int32)
    found = set()
    matched_bg = []
    compfeat = fnt.partial(compute_seq_features, limits['len_bound_hi'])
    zero_cand = 0
    while not timeout.is_set() and len(found) < limits['threshold']:
        # select all free start coordinates
        rand_starts = coordinates[~covered]
        num_candidates = max(25, int(np.ceil((limits['num_samples'] - len(found)) * 1.1)))
        rand_starts = rng.choice(rand_starts, size=num_candidates, replace=False)
        rand_lengths = rng.randint(low=limits['len_bound_lo'],
                                   high=limits['len_bound_hi'],
                                   size=rand_starts.size)
        cand_regions = []
        cand_points = []
        for s, e in zip(rand_starts, rand_starts + rand_lengths):
            if e > chrom_bound or covered[s:e].any():
                continue
            cand_seq = chromseq[s:e]
            cand_feat = compfeat(cand_seq)
            cand_regions.append((s, e, cand_feat[0], cand_feat[1]))
            cand_points.append(cand_feat)
        if len(cand_points) == 0:
            # this can happen, e.g., for small
            # chromosomes like chrM
            chrom_len = float(covered.size)
            if covered.sum() / chrom_len > 0.9:
                # most of the chromosome covered, unlikely
                # to find more matches
                break
            if zero_cand > 20:
                # I assume this is tantamount to the above
                # but just as a safeguard against wasting
                # more CPU cycles
                break
            zero_cand += 1
            continue
        cand_tree = spat.cKDTree(cand_points)
        used_cand = set()
        for eps in np.arange(start=params['epsinit'],
                             stop=params['epsmax'] + params['epsstep'],
                             step=params['epsstep']):
            # bulk query: all candidates against all samples
            neighbors = cand_tree.query_ball_tree(nntree, r=eps, p=np.inf)
            for cand_idx, matches in enumerate(neighbors):
                if cand_idx in used_cand:
                    continue
                sample_idx = [idxmap[m] for m in matches]
                sample_idx = [s for s in sample_idx if s not in found]
                try:
                    found.add(sample_idx[0])
                    used_cand.add(cand_idx)
                    cand_reg = cand_regions[cand_idx]
                    covered[cand_reg[0]:cand_reg[1]] = 1
                    cand_reg = ('', cand_reg[0], cand_reg[1], sample_idx[0], cand_reg[2], cand_reg[3])
                    matched_bg.append(cand_reg)
                except IndexError:
                    continue
        if len(found) >= limits['threshold']:
            break
        nntree, idxmap = remove_matched_regions(nntree, idxmap, found)
    return matched_bg


def remove_matched_regions(nntree, idxmap, found):
    """
    :param nntree:
    :param idxmap:
    :param found:
    :return:
    """
    keep_rows = []
    for k, v in idxmap.items():
        if v not in found:
            keep_rows.append((k, v))
    keep_rows = sorted(keep_rows)
    new_data = nntree.data[[t[0] for t in keep_rows], :]
    new_idxmap = dict([(i, v) for i, (r, v) in enumerate(keep_rows)])
    nntree = spat.cKDTree(new_data)
    return nntree, new_idxmap


def find_background_regions(arguments):
    """
    :param arguments:
    :return:
    """
    chrom, params, fg_regions = arguments
    data_limits = compute_data_limits(fg_regions, params)
    chrom_seq = load_chromosome_sequence(chrom, params['genome'])

    fg_regions, taken, idxmap = bulk_add_seq_features(fg_regions, chrom_seq, data_limits['len_bound_hi'])
    assert fg_regions.shape[0] == data_limits['num_samples'], \
        'Lost foreground entries: {} vs {}'.format(data_limits['num_samples'],
                                                   fg_regions.shape[0])
    kdtree = spat.cKDTree(fg_regions)
    try:
        matched = start_relaxed_search(kdtree, chrom_seq, taken, idxmap, data_limits, params)
    except Exception:
        trb.print_exc(file=sys.stderr)
        return 'Error processing chromosome {}'.format(chrom), []
    else:
        return chrom, matched


def run_background_match(args):
    """
    :param args:
    :return:
    """
    foreground = read_input_file(args.inputfile)
    params = vars(args)
    foreground = [(chrom, params, regions) for chrom, regions in foreground.items()]

    pool = mp.Pool(args.workers)
    try:
        resit = pool.imap_unordered(find_background_regions, foreground)
        out_mode = 'w'
        for chrom, matches in resit:
            if chrom.startswith('Error'):
                sys.stderr.write(chrom)
                continue  # maybe we are lucky and the rest finishes
            if len(matches) == 0:
                sys.stderr.write('\nWarning: no matches for chromosome {}\n'.format(chrom))
                continue
            if args.outputfile == 'stdout':
                [sys.stdout.write(chrom + '\t'.join(map(str, reg)) + '\n') for reg in matches]
            else:
                with open(args.outputfile, out_mode) as dump:
                    [dump.write(chrom + '\t'.join(map(str, reg)) + '\n') for reg in matches]
                out_mode = 'a'
            print('Results from chrom {} saved'.format(chrom))
        pool.close()
        pool.join()
    except:
        pool.close()
        pool.terminate()
        raise
    return 0


if __name__ == '__main__':
    try:
        args = parse_command_line()
        _ = run_background_match(args)
    except Exception as e:
        sys.stderr.write('\nError: {}\n'.format(e))
        trb.print_exc(file=sys.stderr)
        sys.exit(1)
    else:
        sys.exit(0)
