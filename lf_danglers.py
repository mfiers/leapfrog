#!/usr/bin/env python

import sys
import itertools
import subprocess as sp
import argparse


def run_bowtie2(db, fq, preset, threads):
    cl = "bowtie2 --%s --reorder --no-head --quiet -p %d -x %s -U %s" % (
        preset, threads, db, fq)
    process = sp.Popen(cl, shell=True, stdout=sp.PIPE)
    for line in process.stdout:
        if line[0] == '@':
            continue
        yield line.split()


def identify_danglers(arguments):
    i = 0
    j = 1
    with open(arguments.output, 'w') as F:
        for forward, reverse in itertools.izip(run_bowtie2(arguments.bowtie2_database,
                                                           arguments.forward_reads,
                                                           arguments.bowtie_preset,
                                                           arguments.threads),
                                               run_bowtie2(arguments.bowtie2_database,
                                                           arguments.reverse_reads,
                                                           arguments.bowtie_preset,
                                                           arguments.threads)):
            i += 1
            fflag = int(forward[1])
            rflag = int(reverse[1])
            if fflag & 0x4 == rflag & 0x4:
                # either both map, or both do not map
                # not dangling
                continue

            if fflag & 0x4:
                # forward is dangling
                tag = 'R' + {0: '+', 0x10: '-'}[rflag & 0x10]
                F.write('@%s__%s__%s\n' % (reverse[2], tag, forward[0]))
                F.write('%s\n+\n%s\n' % (forward[9], forward[10]))

            else:
                # reverse is dangling
                tag = 'F' + {0: '+', 0x10: '-'}[fflag & 0x10]
                F.write('@%s__%s__%s\n' % (forward[2], tag, reverse[0]))
                F.write('%s\n+\n%s\n' % (reverse[9], reverse[10]))


def parse_args(arguments):
    parser = argparse.ArgumentParser(description='find danglers')
    parser.add_argument('-t', '--threads', default=10, type=int)
    parser.add_argument('bowtie2_database')
    parser.add_argument('forward_reads')
    parser.add_argument('reverse_reads')
    parser.add_argument('output')
    parser.add_argument('-p', '--bowtie_preset', default='fast')
    return parser.parse_args()


def main():
    arguments = parse_args(sys.argv[1:])
    identify_danglers(arguments)


if __name__ == '__main__':
    main()
