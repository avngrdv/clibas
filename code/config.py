experiment = 'unnamed_experiment'

class constants:
    '''
    Star symbol (*) is reserved for stop codons.
    
    Plus and underscore symbols (+ and _) are internally reserved tokens.
    Numerals (1234567890) are internally reserved for library design 
    specification. These symbols (123456790+_) should not be used to
    encode amino acids.
    
    Other symbols are OK.
    '''
    
    translation_table = {
                            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
                       }


    bases = ('T', 'C', 'A', 'G')
    
    #SMILES strings corresponding to all encoded amino acids;
    #sorted alphabetically by amino acid one-letter codes
    aa_SMILES = (   
                 'N[C@@H](C)C(=O)',            
                 'N[C@@H](CS)C(=O)',        
                 'N[C@@H](CC(=O)O)C(=O)',
                 'N[C@@H](CCC(=O)O)C(=O)',
                 'N[C@@H](Cc1ccccc1)C(=O)',
                 'NCC(=O)',
                 'N[C@@H](Cc1c[nH]cn1)C(=O)',
                 'N[C@@H]([C@H](CC)C)C(=O)',
                 'N[C@@H](CCCCN)C(=O)',        
                 'N[C@@H](CC(C)C)C(=O)',
                 'N[C@@H](CCSC)C(=O)',
                 'N[C@@H](CC(=O)N)C(=O)',
                 'O=C[C@@H]1CCCN1',  
                 'N[C@@H](CCC(=O)N)C(=O)',
                 'N[C@@H](CCCNC(=N)N)C(=O)',
                 'N[C@@H](CO)C(=O)',
                 'N[C@@H]([C@H](O)C)C(=O)',
                 'N[C@@H](C(C)C)C(=O)',
                 'N[C@@H](Cc1c[nH]c2c1cccc2)C(=O)',
                 'N[C@@H](Cc1ccc(O)cc1)C(=O)'
                )
    
    #nucleotide complement table;
    #bases are represented by their ascii symbol numbers
    complement_table = {65: 84,
                        67: 71,
                        84: 65, 
                        71: 67, 
                        78: 78}
          
class LibaryDesigns:
        
    #DNA library design parameters     
    dna_templates = [
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                     'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG112112112112112112112112112112112112112112112TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA',
                    ]

    dna_monomers = {
                    1: ('A', 'G', 'T', 'C'),
                    2: ('G', 'T'),
                   }
                    
    #peptide library design parameters
    pep_templates = [
                     '211113GSGSGS',
                     '2111113GSGSGS',
                     '21111113GSGSGS',
                     '211111113GSGSGS',
                     '2111111113GSGSGS',
                     '21111111113GSGSGS',
                     '211111111113GSGSGS',
                     '2111111111113GSGSGS',
                     '21111111111113GSGSGS',
                     '211111111111113GSGSGS',
                     '2111111111111113GSGSGS',
                     '21111111111111113GSGSGS',
                    ]
            
    pep_monomers = {
                    1: ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
                        'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'),
                    2: ('M'),
                    3: ('C')
                   }
    
class TrackerConfig:
    
    #directory holding sequencing data files (fastq or fastq.gz)
    seq_data = '../sequencing_data'
        
    #directory for writing logs to
    logs = '../logs'
    
    #directory that stores fastqparser outputs
    parser_out = '../parser_outputs'
    
    #directory that stores outputs of data analysis ops
    analysis_out = '../data_analysis'

class LoggerConfig:
       
    #verbose loggers print to the console
    verbose = True
    
    #write logs to file
    log_to_file = True

    #logger level; accepted any of ('INFO', 'WARNING', 'ERROR')
    level = 'INFO'    
    
class FastqParserConfig:
    
    #a RE pattern that has to match in order to initiate the orf
    #it will be used when translation is performed 
    #with force_translation=False flag
    utr5_seq = 'AGGAGAT......ATG'