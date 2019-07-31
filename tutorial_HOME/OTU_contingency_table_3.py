#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read all fasta files and build a sorted OTU contingency
    table. Usage: python OTU_contingency_table.py [input files]
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr> modified by Benoit Perez-Lamarque <benoit.perez@ens.fr>"
__date__ = "2019/03/07"
__version__ = "$Revision: 5.0"

import os
import re
import sys
import operator

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def representatives_parse():
    """
    Get seed sequences.
    """
    representatives_file = sys.argv[1]
    representatives = dict()
    with open(representatives_file, "rU") as representatives_file:
        for line in representatives_file:
            if line.startswith(">"):
                amplicon = line.strip(">;\n")
            else:
                representatives[amplicon] = line.strip("\n")

    return representatives


def stats_parse_old():
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = sys.argv[2]
    stats = dict()
    seeds = dict()
    with open(stats_file, "rU") as stats_file:
        for line in stats_file:
            cloud, mass, seed, seed_abundance = line.strip().split(separator)[0:4]
            stats[seed] = int(mass)
            seeds[seed] = (int(seed_abundance), int(cloud))
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.items(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds



def stats_parse():
    """
    Map OTU seeds and stats.
    """
    separator = "   "
    stats_file = sys.argv[2]
    stats = dict()
    seeds = dict()
    with open(stats_file, "rU") as stats_file:
        for line in stats_file:
            cloud, seed = line.strip().split(separator)[0:2]
            stats[cloud] = 1
            seeds[cloud] = (1, 1)
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.items(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds


def swarms_parse():
    """
    Map OTUs.
    """
    separator = "\t"
    swarms_file = sys.argv[3]
    swarms = dict()
    with open(swarms_file, "rU") as swarms_file:
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::]
            seed = amplicons[0]
            swarms[seed] = [amplicons[1::]]

    return swarms


def uchime_parse():
    """
    Map OTU's chimera status.
    """
    separator = " "
    uchime_file = sys.argv[4]
    uchime = dict()
    with open(uchime_file, "rU") as uchime_file:
        for line in uchime_file:
            OTU = line.strip().split("\t")
            try:
                seed = OTU[1].split(";")[0]
            except IndexError:  # deal with partial line (missing seed)
                continue
            try:
                status = OTU[17]
            except IndexError:  # deal with unfinished chimera detection runs
                status = "NA"
            uchime[seed] = status

    return uchime


def stampa_parse():
    """
    Map amplicon ids and taxonomic assignments.
    """
    separator = "\t"
    stampa_file = sys.argv[5]
    stampa = dict()
    with open(stampa_file, "rU") as stampa_file:
        for line in stampa_file:
            amplicon, identity, taxonomy = line.strip().split(separator)
            stampa[amplicon] = (identity, taxonomy)

    return stampa


def fasta_parse():
    """
    Map amplicon ids, abundances and samples.
    """
    fasta_files = sys.argv[6:]
    samples = dict()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        sample = os.path.basename(fasta_file)
        sample = os.path.splitext(sample)[0]
        samples[sample] = samples.get(sample, 0) + 1
        with open(fasta_file, "rU") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon = line.strip(">;\n")
                    abundance = 1
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        # deal with duplicated samples
                        amplicons2samples[amplicon][sample] = amplicons2samples[amplicon].get(sample, 0) + abundance
    # deal with duplicated samples
    duplicates = [sample for sample in samples if samples[sample] > 1]
    if duplicates:
        print("Warning: some samples are duplicated", file=sys.stderr)
        print("\n".join(duplicates), file=sys.stderr)
    samples = sorted(samples.keys())

    return amplicons2samples, samples


def print_table(representatives, stats, sorted_stats,
                swarms, uchime, amplicons2samples,
                samples, seeds,
                stampa):
    """
    Export results.
    """
    # Print table header
    print("OTU", "total", "cloud",
          "amplicon", "length", "abundance",
          "chimera", "spread",
          "sequence", "identity", "taxonomy",
          "\t".join(samples),
          sep="\t", file=sys.stdout)

    # Print table content
    i = 1
    for seed, abundance in sorted_stats:
        sequence = representatives[seed]
        occurrences = dict([(sample, 0) for sample in samples])
        for amplicons in swarms[seed]:
            for amplicon in amplicons:
                for sample in samples:
                    occurrences[sample] += amplicons2samples[amplicon].get(sample, 0)
        spread = len([occurrences[sample] for sample in samples if occurrences[sample] > 0])
        sequence_abundance, cloud = seeds[seed]

        
        # Chimera checking (deal with incomplete cases. Is it useful?)
        if seed in uchime:
            chimera_status = uchime[seed]
        else:
            chimera_status = "NA"

        # Chimera checking (deal with incomplete cases. Is it useful?)
        if seed in stampa:
            identity, taxonomy = stampa[seed]
        else:
            identity, taxonomy = "NA", "NA"

        # output
        print(i, abundance, cloud,
              seed, len(sequence), sequence_abundance,
              chimera_status, spread, sequence,
              identity, taxonomy,
              "\t".join([str(occurrences[sample]) for sample in samples]),
              sep="\t", file=sys.stdout)
        i += 1

    return


def main():
    """
    Read all fasta files and build a sorted OTU contingency table.
    """
    #print("Parse taxonomic results")
    representatives = representatives_parse()

    #print("Parse stats")
    stats, sorted_stats, seeds = stats_parse()

#print("Parse swarms")
    swarms = swarms_parse()

#print("Parse uchime")
    uchime = uchime_parse()

    # Parse taxonomic assignment results
    stampa = stampa_parse()
    
    #print("Parse fasta files")
    amplicons2samples, samples = fasta_parse()

#print("Print table header")
    print_table(representatives, stats, sorted_stats,
            swarms,uchime, amplicons2samples, samples,
            seeds,
            stampa)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
