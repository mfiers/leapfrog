#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
import collections

def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('-a', '--min_coverage', type=int, default=4,
                        help=("Minimal coverage for a family to be "
                              "exported as a peak - the number "
                              "checked is the maximal coverage for "
                              "that set of mapping reads. (default 4)"))
    parser.add_argument('-m', '--min_diff', type=int, default=5,
                        help=("Stringency to determine if a hit maps "
                              "uniquely. This score is the "
                              "minimal allowed difference between the "
                              "scores of the first and second "
                              "hit before a read is assumed to be "
                              "mapped to it's correct, unique, "
                              "location. (default=5)"))
    parser.add_argument('-c', '--trim_cov', type=float, default=3.,
                        help=("Minimal coverage of a read region "
                              "that is considered part of the "
                              "peak. This number is used to trim "
                              "trailing edges from a peak, or, if "
                              "necessary to cut a peak into an "
                              "number of smaller peaks. Note, if you "
                              "set a value below 0, it will be "
                              "interpreted as a fraction of the "
                              "maximal coverage observed for that "
                              "peak. (default: 3)"))
    parser.add_argument('-C', '--min_trim_cov', type=int, default=3,
                        help=("When using a fraction for -c, this "
                              "value sets a lower limit - anything "
                              "below this value will be trimmed "
                              "for sure. (default 3)"))
    parser.add_argument('-q', '--output_nonunique',
                        action='store_true', default=False,
                        help=("also output peaks which are no likely "
                              "to be uniquely mapped"))
    return parser.parse_args(args)


def extract_clusters(sam):
    """
    Extracts all clusters of overlapping reads from a sorted, indexed pysam object.
    Generates a dictionary per cluster.
    Accepts a pysam object.
    Returns a dictionary generator.
    """

    cluster = {"reads": [],
               "chromosome": "",
               "start": -1,
               "stop": -1,}

    for i, read in enumerate(sam.fetch()):
        read_chromosome = sam.getrname(read.tid)
        read_start = read.pos
        read_stop = read.pos + read.qlen

        # if read overlaps the current cluster
        if (read_chromosome == cluster["chromosome"]) and (read_start < cluster["stop"]):

            # add the read to the current cluster
            cluster["reads"].append(read)
            cluster["stop"] = read_stop

        # else read is the start of a new cluster
        else:

            # yield the previous cluster but skip the first blank cluster
            if i > 0:
                yield cluster

            # create a new cluster dictionary based on the current read
            cluster["reads"] = [read]
            cluster["chromosome"] = read_chromosome
            cluster["start"] = read_start
            cluster["stop"] = read_stop

    # ensure the final cluster is not skipped
    yield cluster


def split_families(cluster_generator):
    """
    Subdivides read-clusters based on read family.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for cluster in cluster_generator:
        families = collections.defaultdict(list)
        for read in cluster["reads"]:
            try:
                family = read.qname.split('__')[0].rsplit('/', 1)[1]
            except IndexError:
                family = read.qname.split('__')[0]
            families[family].append(read)

        for family, reads in families:
            start = min([read.pos for read in reads])
            stop = max([read.pos for read in reads])
            sub_cluster = {"reads": reads,
                           "chromosome": cluster["chromosome"],
                           "family": family,
                           "start": start,
                           "stop": stop}
            yield sub_cluster


def split_orientation():
    """
    Subdivides read-clusters based on read orientation.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def map_clusters():
    """
    Appends a numpy array of read depth.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def trim_clusters():
    """
    Trims read-clusters based on a read depth.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def split_gaps():
    """
    Subdivides read-clusters with a read depth of 0.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def filter_unique():
    """
    Filters read-clusters by the amount of uniquely mapped reads.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def filter_depth():
    """
    Filters read-clusters based on maximum read depth.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    pass


def construct_gff():
    """
    Accepts a dictionary generator.
    Prints a gff file to standard out.
    """
    pass


def main():
    pass


if __name__ == '__main__':
    main()