#!/usr/bin/env python

#################################
# Copied from "Dina" at https://github.com/dinovski
# on gist https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
# Dina's code looked like the best way to computer N50 (which is surprisingly not in a library
# anywhere).  I've modified the code to serve as a module inside DDV.
#
# The code was further modified to remove a numpy dependency and keep it lightweight.

# calculate N50 from fasta file
# N50 = contig length so that half of the contigs are longer and 1/2 of contigs are shorter
from __future__ import print_function, division, absolute_import, with_statement
import sys
from itertools import chain

from DNASkittleUtils.Contigs import read_contigs


def cumulative_sum(numbers_list):
    running_sums = []
    current_sum = 0
    for i in numbers_list:
        current_sum += i
        running_sums.append(int(current_sum))
    return running_sums


def collect_n50_stats(scaffold_lengths, prefix=''):
    """N50:
    the length of the shortest contig such that the sum of contigs of equal
    length or longer is at least 50% of the total length of all contigs"""

    stats = {}

    # sort contigs longest>shortest
    all_len = sorted(scaffold_lengths, reverse=True)
    csum = cumulative_sum(all_len)

    assembly_size = sum(scaffold_lengths)
    stats[prefix + 'N'] = int(assembly_size)
    halfway_point = (assembly_size // 2)

    # get index for cumsum >= N/2
    for i, x in enumerate(csum):
        if x >= halfway_point:
            stats[prefix + 'N50'] = all_len[i]
            break

    # N90
    stats[prefix + 'nx90'] = int(assembly_size * 0.90)

    # index for csumsum >= 0.9*N
    for i, x in enumerate(csum):
        if x >= stats[prefix + 'nx90']:
            stats[prefix + 'N90'] = all_len[i]
            break

    return stats


def scaffold_lengths_from_fasta(input_fasta_path):
    scaffolds = read_contigs(input_fasta_path)
    lengths = [len(x.seq) for x in scaffolds]
    return scaffolds, lengths


def split_by_N(scaffolds):
    length_collection = set()
    for scaffold in scaffolds:
        pieces = scaffold.seq.split('N')
        length_collection.add((len(p) for p in pieces))
    lengths = list(chain(*length_collection))
    return lengths


def all_stats(input_fasta):
    scaffolds, lengths = scaffold_lengths_from_fasta(input_fasta)
    scaffold_stats = collect_n50_stats(lengths, prefix='Scaffold ')
    contig_lengths = split_by_N(scaffolds)
    contig_stats = collect_n50_stats(contig_lengths, prefix='Contig ')
    scaffold_stats.update(contig_stats)
    scaffold_stats['N%'] = (1 - (scaffold_stats['Contig N'] / float(scaffold_stats['Scaffold N']))) * 100
    return scaffold_stats


if __name__ == '__main__':
    input_fasta_name = sys.argv[1]
    assembly_stats = all_stats(input_fasta_name)
    label_order = ['Scaffold N', 'Scaffold N50', 'Scaffold N90', 'Scaffold nx90',
                   'Contig N', 'Contig N50', 'Contig N90', 'Contig nx90',
                   'N%']
    for key in label_order:
        print(key + ":", "{:,}".format(assembly_stats[key]))
    for key in assembly_stats:  # unordered labels
        if key not in label_order:
            print(key + ":", "{:,}".format(assembly_stats[key]))


