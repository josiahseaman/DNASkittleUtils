from __future__ import print_function, division, absolute_import, with_statement

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

