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


def __do_write(filestream, seq, header=None):
    """Specialized function for writing sets of headers and sequence in FASTA.
    It chunks the file up into 70 character lines, but leaves headers alone"""
    if header is not None:
        filestream.write(header + '\n')  # double check newlines
    try:
        for line in chunks(seq, 70):
            filestream.write(line + '\n')
    except Exception as e:
        print(e)


def _write_fasta_lines(filestream, seq):
    import _io
    assert isinstance(filestream,
                      _io.TextIOWrapper)  # I'm actually given a file name and have to open it myself
    contigs = seq.split('\n')
    index = 0
    while index < len(contigs):
        if len(contigs) > index + 1 and contigs[index].startswith('>') and contigs[index + 1].startswith('>'):
            print("Warning: Orphaned header:", contigs[index])
        if contigs[index].startswith('>'):
            header, contents = contigs[index], contigs[index + 1]
            index += 2
        else:
            header, contents = None, contigs[index]
            index += 1
        __do_write(filestream, contents, header)


def write_complete_fasta(file_path, seq_content_array, header=None):
    """This function ensures that all FASTA files start with a >header\n line"""
    with open(file_path, 'w') as filestream:
        if seq_content_array[0] != '>':  # start with a header
            temp_content = seq_content_array
            if header is None:
                header = '>%s\n' % just_the_name(file_path)
            if isinstance(temp_content, list):
                seq_content_array = [header]
            else:
                seq_content_array = array('u', header)
            seq_content_array.extend(temp_content)
        _write_fasta_lines(filestream, ''.join(seq_content_array))


def write_contigs_to_file(out_filename, contigs):
    with open(out_filename, 'w') as outfile:
        for contig in contigs:
            __do_write(outfile, header='>' + contig.name, seq=contig.seq)
    print("Done writing ", len(contigs), "contigs and {:,}bp".format(sum([len(x.seq) for x in contigs])))


class BlankIterator:
    def __init__(self, filler):
        self.filler = filler

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.filler * (index.stop - index.start)
        else:
            return self.filler
