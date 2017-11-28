from __future__ import print_function, division, absolute_import, with_statement

import os
import re as regex
import shutil
import textwrap
from array import array
from collections import namedtuple

from DNASkittleUtils.CommandLineUtils import just_the_name

Batch = namedtuple('Batch', ['chr', 'fastas', 'output_folder'])
nucleotide_complements = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N', 'X': 'X'}


def complement(plus_strand):
    return nucleotide_complements[plus_strand]


def rev_comp(plus_strand):
    return ''.join([nucleotide_complements[a] for a in reversed(plus_strand)])


class ReverseComplement:
    def __init__(self, seq, annotation=False):
        """Lazy generator for being able to pull out small reverse complement
        sections out of large chromosomes"""
        self.seq = seq
        self.length = len(seq)
        self.annotation = annotation

    def __getitem__(self, key):
        if isinstance(key, slice):
            end = self.length - key.start
            begin = self.length - key.stop
            if end < 0 or begin < 0 or end > self.length:
                raise IndexError("%i %i vs. length %i" % (end, begin, self.length))
            piece = self.seq[begin: end]
            return rev_comp(piece) if not self.annotation else ''.join(reversed(piece))
        letter = self.seq[self.length - key - 1]
        return complement(letter) if not self.annotation else letter

    def __len__(self):
        return 0


def pretty_contig_name(contig, title_width, title_lines):
    """Since textwrap.wrap break on whitespace, it's important to make sure there's whitespace
    where there should be.  Contig names don't tend to be pretty."""
    pretty_name = contig.name.replace('_', ' ').replace('|', ' ').replace('chromosome chromosome',
                                                                          'chromosome')
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)  # don't ask
    if title_width < 20 and len(
            pretty_name) > title_width * 1.5:  # this is a suboptimal special case to try and
        # cram more characters onto the two lines of the smallest contig titles when there's not enough space
        # For small spaces, cram every last bit into the line labels, there's not much room
        pretty_name = pretty_name[:title_width] + '\n' + pretty_name[title_width:title_width * 2]
    else:  # this is the only case that correctly bottom justifies one line titles
        pretty_name = '\n'.join(textwrap.wrap(pretty_name, title_width)[:title_lines])  # approximate width
    return pretty_name


def copytree(src, dst, symlinks=False, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst, exist_ok=True)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)


def chunks(seq, size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def first_word(string):
    import re
    if '\\' in string:
        string = string[string.rindex('\\') + 1:]
    return re.split('[\W_]+', string)[0]


class BlankIterator:
    def __init__(self, filler):
        self.filler = filler

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.filler * (index.stop - index.start)
        else:
            return self.filler
