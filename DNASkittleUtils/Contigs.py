from __future__ import print_function, division, absolute_import, with_statement
from array import array

from DNASkittleUtils.CommandLineUtils import just_the_name
from DNASkittleUtils.DDVUtils import chunks


class Contig:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __repr__(self):
        return '< "%s" %i nucleotides>' % (self.name, len(self.seq))


def read_contigs(input_file_path):
    contigs = []
    current_name = ""
    seq_collection = []

    # Pre-read generates an array of contigs with labels and sequences
    with open(input_file_path, 'r') as streamFASTAFile:
        for read in streamFASTAFile.read().splitlines():
            if read == "":
                continue
            if read[0] == ">":
                # If we have sequence gathered and we run into a second (or more) block
                if len(seq_collection) > 0:
                    sequence = "".join(seq_collection)
                    seq_collection = []  # clear
                    contigs.append(Contig(current_name, sequence))
                current_name = read[1:]  # remove >
            else:
                # collects the sequence to be stored in the contig, constant time performance don't concat strings!
                seq_collection.append(read.upper())

    # add the last contig to the list
    sequence = "".join(seq_collection)
    contigs.append(Contig(current_name, sequence))
    return contigs


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
    assert isinstance(filestream, _io.TextIOWrapper)  # I'm actually given a file name and have to open it myself
    contigs = seq.split('\n')
    index = 0
    while index < len(contigs):
        if len(contigs) > index + 1 and contigs[index].startswith('>') and contigs[index+1].startswith('>'):
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
    if hasattr(contigs, '__len__'):  # in case of iterator e.g. itertools.chain()
        print("Done writing ", len(contigs), "contigs and {:,}bp".format(sum([len(x.seq) for x in contigs])))
    else:
        print("Done writing {:,}bp".format(sum([len(x.seq) for x in contigs])))

def pluck_contig(chromosome_name, genome_source):
    """Scan through a genome fasta file looking for a matching contig name.  When it find it, find_contig collects
    the sequence and returns it as a string with no cruft."""
    chromosome_name = '>' + chromosome_name
    print("Searching for", chromosome_name)
    seq_collection = []
    printing = False
    with open(genome_source, 'r') as genome:
        for line in genome:
            if line.startswith('>'):
                # headers.append(line)
                line = line.rstrip()
                if line.upper() == chromosome_name.upper():
                    printing = True
                    print("Found", line)
                elif printing:
                    break  # we've collected all sequence and reached the beginning of the next contig
            elif printing:  # This MUST come after the check for a '>'
                line = line.rstrip()
                seq_collection.append(line.upper())  # always upper case so equality checks work
    if not len(seq_collection):
        # File contained these contigs:\n" + '\n'.join(headers)
        raise IOError("Contig not found." + chromosome_name + "   inside " + genome_source)
    return ''.join(seq_collection)

