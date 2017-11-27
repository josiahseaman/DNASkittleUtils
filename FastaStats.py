#!/usr/bin/env python

#################################
# Copied from "Dina" at https://github.com/dinovski
# on gist https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
# Dina's code looked like the best way to computer N50 (which is surprisingly not in a library
# anywhere).  I've modified the code to serve as a module inside DDV.

# calculate N50 from fasta file
# N50 = contig length so that half of the contigs are longer and 1/2 of contigs are shorter
from __future__ import print_function, division
import sys

from Contigs import read_contigs


def cumulative_sum(numbers_list):
    running_sums = []
    current_sum = 0
    for i in numbers_list:
        current_sum += i
        running_sums.append(int(current_sum))
    return running_sums


def collect_n50_stats(scaffold_lengths):
    """N50:
    the length of the shortest contig such that the sum of contigs of equal
    length or longer is at least 50% of the total length of all contigs"""

    stats = {}

    # sort contigs longest>shortest
    all_len = sorted(scaffold_lengths, reverse=True)
    csum = cumulative_sum(all_len)

    assembly_size = sum(scaffold_lengths)
    print("N: %d" % int(assembly_size))
    halfway_point = (assembly_size // 2)

    # get index for cumsum >= N/2
    for i, x in enumerate(csum):
        if x >= halfway_point:
            stats['N50'] = all_len[i]
            print("N50: %s" % stats['N50'])
            break

    # N90
    stats['nx90'] = int(assembly_size * 0.90)

    # index for csumsum >= 0.9*N
    for i, x in enumerate(csum):
        if x >= stats['nx90']:
            stats['N90'] = all_len[i]
            print("N90: %s" % stats['N90'])
            break

    return stats


def scaffold_lengths_from_fasta(input_fasta_path):
    contigs = read_contigs(input_fasta_path)
    lengths = [len(x.seq) for x in contigs]
    return lengths


if __name__ == '__main__':
    input_fasta_name= sys.argv[1]
    lengths = scaffold_lengths_from_fasta(input_fasta_name)
    assembly_stats = collect_n50_stats(lengths)
    for key in assembly_stats:
        print(key, "{:,}".format(assembly_stats[key]))
