import os
from utils.lib_design import LibraryDesign

experiment = 'unnamed_experiment'


class constants:
    '''
    Underscore symbol (_) is reserved for stop codons;
    Plus symbol (+) is an internally reserved token;
    Numerals (1234567890) are internally reserved for
    library design specification. Other symbols can
    be used for encoding.

    Although the class holds multiple attributes, only
    the codon table should be edited. Everything else
    is inferred automatically from it.
    '''
    # NNU15MGS3
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'B', 'TGT': 'B', 'TGA': '*', 'TGG': 'W',
    }

    bases = ('T', 'C', 'A', 'G')
    complement_table = str.maketrans('ACTGN', 'TGACN')

    global _reserved_aa_names
    _reserved_aa_names = ('_', '+', '1', '2', '3', '4',
                          '5', '6', '7', '8', '9', '0')

    aas = tuple(
        sorted(set(x for x in codon_table.values() if x not in _reserved_aa_names)))
    codons = tuple(sorted(set(x for x in codon_table.keys())))
    aa_dict = {aa: i for i, aa in enumerate(aas)}


class ParserConfig:

    # a RE pattern that has to match in order to initiate the orf
    # is only used when force_translation == False
    utr5_seq = 'AGGAGAT......ATG'

    # DNA library design
    D_design = LibraryDesign(

        # NNU15MGS3
        templates=[
                  'TAATACGACTCACTATAGGGTTGAACTTTAAGTAGGAGATATATCCATG112112112112112112112112112112112112112112112ATGGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
        ],

        monomers={
            1: ('A', 'G', 'T', 'C'),
            2: ('T'),
        },

        lib_type='dna'

    )

    # peptide library design
    P_design = LibraryDesign(

        # NNU15MGS3
        templates=[
                   'M111111111111111MGSGSGS',
        ],

        monomers={
            1: ('A', 'B', 'D', 'F', 'G', 'H', 'I', 'L', 'N', 'P', 'R', 'S', 'T', 'V', 'Y'),
        },

        lib_type='pep'

    )


class TrackerConfig:

    # directory holding sequencing data files (fastq or fastq.gz)
    seq_data = '../sequencing_data'

    # directory for writing logs to
    logs = '../logs'

    # directory that stores fastqparser outputs
    parser_out = '../parser_outputs'


class LoggerConfig:

    # logger name
    name = experiment + '_logger'

    # verbose loggers print to the console
    verbose = True

    # write logs to file
    log_to_file = True

    # log filename; when None, the name will be inferred by the Logger itself
    log_fname = os.path.join(TrackerConfig.logs, experiment + ' logs')
