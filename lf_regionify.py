#!/usr/bin/env python

import numpy as np
import pysam
import logging
import argparse
import collections

logging.basicConfig()
lg = logging.getLogger(__name__)

md_help = " ".join("""
Stringency to determine if a hit maps uniquely. This score is the
minimal allowed difference between the scores of the first and second
hit before a read is assumed to be mapped to it's correct, unique,
location. (default=5)
""".split())

mincov_help = " ".join("""
Minimal coverage for a family to be exported as a peak - the number
checked is the maximal coverage for that set of mapping
reads. (default 4)
""".split())

peakfilter_help = " ".join("""
Minimal coverage of a read region that is considered part of the
peak. This number is used to trim trailing edges from a peak, or, if
necessary to cut a peak into an number of smaller peaks. Note, if you
set a value below 0, it will be interpreted as a fraction of the
maximal coverage observed for that peak. (default: 3)
""".split())


mintrimcov_help = " ".join("""
When using a fraction for -c, this value sets a lower limit - anything
below this value will be trimmed for sure. (default 3)
""").split()

parser = argparse.ArgumentParser('Identify transposon flanking regions')

parser.add_argument('input_bam')

parser.add_argument('-a',
                    '--min_coverage',
                    type=int, default=4,
                    help=mincov_help)

parser.add_argument('-m',
                    '--min_diff',
                    type=int,
                    default=5,
                    help=md_help)

parser.add_argument('-c',
                    '--trim_cov',
                    type=float,
                    default=3.,
                    help=peakfilter_help)

parser.add_argument('-C',
                    '--min_trim_cov',
                    type=int,
                    default=3,
                    help=mintrimcov_help)

parser.add_argument('-q',
                    '--output_nonunique',
                    action='store_true',
                    default=False,
                    help='also output peaks which are no likely to be uniquely mapped')

args = parser.parse_args()

sam = pysam.Samfile(args.input_bam, 'rb')

bump_chrom = ""
bump_start = -1
bump_stop = -1
bump_reads = []

FC = 0


def famdump(chrom, family, reads):
    global FC
    FC += 1

    no_unique_hits = 0
    for r in reads:
        AS = -999
        XS = -999
        SDF = 0
        for t in r.tags:
            if t[0] == 'AS':
                AS = t[1]
            if t[0] == 'XS':
                XS = t[1]

        SDF = AS - XS
        # lg.critical("%s %s %s" % ( AS, XS, SDF))
        if SDF >= args.min_diff:
            no_unique_hits += 1

    start = min([r.pos for r in reads])
    stop = max([r.pos + r.qlen for r in reads])

    # generate a regional coverage plot
    covarr = np.zeros((stop-start))
    for r in reads:
        covarr[(r.pos - start):(r.pos + r.qlen - start)] += 1

    # max coverage
    maxcov = max(covarr)

    # set anything below args.trim_cov to 0
    trim_cutoff = args.trim_cov
    if trim_cutoff < 1:
        trim_cutoff *= maxcov
        trim_cutoff = max(args.min_trim_cov, trim_cutoff)

    covarr[covarr < trim_cutoff] = 0

    # trim the peak - start
    while len(covarr) > 0 and covarr[0] == 0:
        start += 1
        covarr = covarr[1:]
    # trim the peak -end
    while len(covarr) > 0 and covarr[-1] == 0:
        stop -= 1
        covarr = covarr[:-1]

    if len(covarr) == 0:
        return 0

    avg_coverage = sum([r.qlen for r in reads]) / float(stop - start)
    if maxcov < args.min_coverage:
        return 0

    if 0 in covarr:
        # go into gapsplitting mode
        gapsplit_rv = 0

        # print "SPLIT"
        # print 'start, stop, maxcov', start, stop, maxcov
        # print covarr
        # after trimming start & stop we should always start inside a
        # peak
        def _process_peak(peakstart, peakstop):
            """
            process a peak - (iteratively send it of for processing

            """
            # identify reads that overlap with this peak
            peakreads = []
            for r in reads:
                # print r.pos, r.pos + r.qlen, start + peakstart, stop + peakstop
                if r.pos > start + peakstop:
                    continue
                if r.pos + r.qlen < start + peakstart:
                    continue
                peakreads.append(r)
            # print 'peak %d to %d with %d reads out of %d ' % (
            # peakstart, peakstop, len(peakreads), len(reads))
            return famdump(chrom, family, peakreads)

        inpeak, peakstart, peakstop = True, 0, 0
        for i, ingap in enumerate(covarr == 0):
            if not ingap:
                # we're not in a gap (yet)
                peakstop = i
                if not inpeak:
                    # and we'er just coming into a peak
                    inpeak = True
                    peakstart = i
            else:
                # we are in a gap -
                if inpeak:
                    # we just left a peak
                    gapsplit_rv += _process_peak(peakstart, peakstop)
                    covarr[peakstart:peakstop+1] = 99
                inpeak = False
        # we're leaving a peak for sure now
        if inpeak:
            gapsplit_rv += _process_peak(peakstart, peakstop)
            covarr[peakstart:peakstop+1] = 99

        return gapsplit_rv

    names = [r.qname.split('__')[0].rsplit('/', 1)[0] for r in reads]
    namecount = reversed(sorted([(float(names.count(e)) / len(names), e) for e in set(names)]))
    namecount = [x for x in namecount if x[0] > 0.1]
    names = ",".join([n[1] for n in namecount])

    fracreverse = float(len([1 for r in reads if r.is_reverse])) / len(reads)
    if fracreverse < 0.1:
        strand = '+'
    elif fracreverse < 0.9:
        strand = '.'
    else:
        strand = '-'
    visstrand = ''.join(sorted([{True:'-', False:'+'}[r.is_reverse] for r in reads]))

    if no_unique_hits >= args.min_diff:
        state = 'UNIQUE'
    else:
        state = 'NOTUNIQ'

    if (not args.output_nonunique) and (state == 'NOTUNIQ'):
        return 0

    print '%s\tREFS\tREFS.%s.%s\t%d\t%d\t%.2f\t%s\t.\tID=reps_%s_%s_%s;Name="%s"' % (
        chrom, state, family, start, stop, avg_coverage, strand, chrom, start, family,
        names)

    return 1


def clean_family_name(f):
    return f.replace('(', '').replace(')', '').replace('#', '')


def bumpdump(chrom, reads):

    fams = collections.defaultdict(list)

    # first subdivide into families
    for read in reads:
        try:
            family = read.qname.split('__')[0].rsplit('/', 1)[1]
        except IndexError:
            family = read.qname.split('__')[0]

        fams[family].append(read)

    rv = 0
    for f in fams:
        rv += famdump(chrom, f, fams[f])
    return rv

bump_count = 0
for i, read in enumerate(sam.fetch()):

    chrom = sam.getrname(read.tid)
    start = read.pos
    stop = read.pos + read.qlen

    if bump_chrom != chrom or start > bump_stop:
        # new bump - process old bump

        if bump_reads:
            bump_count += bumpdump(bump_chrom, bump_reads)
        bump_chrom = chrom
        bump_start = start
        bump_stop = stop
        bump_reads = [read]

    else:
        # still in bump
        # we demand that the reads are sorted - so only need to check the

        bump_stop = max(stop, bump_stop)
        bump_reads.append(read)
