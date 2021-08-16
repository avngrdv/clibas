# FastqProcessor
### An NGS data processor tailored for mRNA display selection workflows

FastqProcessor is python package to analyze NGS sequencing data. The tool is specifically tailored toward analysis of screening outcomes for the selections/screens of genetically encoded peptide/protein libraries (i.e. mRNA , phage or yeast display libraries).\
\
The primary purpose of FastqProcessor is library design matching. Given a particular library design for a genetically encoded peptide/protein library, parse and filter .fastq files to remove noise, low confidence reads, sequences that are too short, too long and much more.\
\
Concomitantly, the tool can collect basic statistics about the data, count sequences, etc.
### Motivation
After working with mRNA-display libraries, I found myself writing and rewriting fastq parser scripts all the time, to accomodate for different workflows. Some libraries encode peptides of a single, fixed length, other don't. Sometimes it is a de novo designed library, and sometimes it's a mutagenesis experiment. I reused my own code, but still had to write entire scripts from scratch for every different application. This package is an attempt to solve this issue by writing a single powerfull parser that can be easily customized for any reasonable application. The goal was to write a package once, and then use it for various applications to create concise (less than 20 lines of code) and readable analysis pipelines. 
# Dependencies
The code was written and tested for:\
\
python 3.8.5\
numpy 1.19.5\
pandas 1.2.4\
matplotlib 3.3.2\
\
The easiest way to get these up and running is downloading and installing [Anaconda individual edition](https://www.anaconda.com/products/individual/).

# Usage
The package provides the tools for analysis. Specific pipelines need to be built for individual applications as explained below. code/config.py and code/main.py also provide examples. _The code is released under the GNU General Public License v3.0 without any explicit or implied warranty. Use it at your own discretion._

# Documentation
## Overall workflow
To analyze a set of .fastq or .fastq.zg files, first a config file (config.py) needs to be edited/created. config.py should hold information about the DNA -> protein translation table (which can be customized), library design information (at both the peptide and DNA level) and other miscellaneous pointers, like the location of input and output folders, logging parameters, etc. \
\
Then, a .py script needs to be written to initialize the parser, queue the operations on data and run the process. The easiest way to initialize the parser and pipeline objects (pipeline is used to enqueue parser operations and run them) is via Dispatcher as shown below:

    from utils.ProcessHandlers import Pipeline, FastqParser
    from utils.Dispatcher import Dispatcher
    
    import config
    
    dispatcher = Dispatcher(config) 
    handlers = (Pipeline, FastqParser)
    
    pip, par = dispatcher.dispatch_handlers(handlers)

FastqParser holds all methods for input/output and data manipulation. These methods can be enqueued in the pipeline. As an example, to enqueue data fetching and DNA-to-peptide translation ops:

    pip.enque([par.fetch_gz_from_dir(), par.translate()])

This will add two ops to the pipeline, but no data fetching or translation will take place just yet. For some ops that contain additional meta parameters, parameter validation will also be performed at this time. Finally, running 

    pip.run(save_summary=True)

will execute the enqueued ops, one by one, in the specified order, passing data from one op to the next. That’s it. If any outputs are expected, they will be written to a directory specified 
in config.py; a log file may also be optionally written if desired.

## LibraryDesign
LibraryDesign is a key object that specifies what kind of library the parser should expect. The object can hold information about arbitrary DNA and peptide libraries using a unified logic as follows. Randomized amino acids/bases (hereafter _tokens_) are indicated as numerals (0-9), whereas tokens which are not subject to randomization (linker sequences, etc) are indicated using the standard one letter encoding (A, C, T, G for DNA), (A, C, D etc for peptides; the encoding must make sense according to the translation table in config.py). A continuous stretch of either random or fixed tokens makes up a “region” in the template sequence. For example:

                seq:      ACDEF11133211AWVFRTQ12345YTPPK
             region:      [-0-][---1--][--2--][-3-][-4-]
        is_variable:      False  True   False True False

Region assignments are made automatically; in the example above, the library contains 5 regions; 3 are “constant regions” and 2 are “variable regions”. Numerals used for variable region tokens should be defined, with one number corresponding to a particular token set. For example, an NNK codon encodes all 20 amino acids, whereas an NNC codon only 15. Thus, all amino acids derived from NNK codons should be encoded by one number, and another number for NNC-encoded positions. LibraryDesign can take several templates of different length to encode libraries with variable regions of variable size. Below is an example of a LibraryDesign initialization:

    lib = LibraryDesign(
        
                    templates=[
                                '211113GSGSGS',
                                '2111113GSGSGS',
                                '21111113GSGSGS',
                              ],
            
                    monomers={
                              1: ('A', 'C', 'D', 'E', 'F', 'G', 'H'),
                              2: ('M'),
                              3: ('C')
                             },
                    
                    lib_type='pep'
                        
                       )

Note that variable positions can encode a single amino acid (amino acids 2 and 3). In this way, there is a considerable flexibility in how a particular library can be represented. When initializing LibraryDesign objects, several rules must be followed:

1.	The topology of every passed template must be identical. Topology is the total number of regions, and the total number of variable regions. Essentially, the templates should only differ in the internal composition of variable regions.
2.	All variable region monomers should be encoded in the translation table (or be one of the four standard DNA bases for DNA libraries; bases N, K etc should be converted to numerals).
3.	Two LibraryDesign objects should be created for the parser (lib_type=’dna’ and lib_type=’pep’). They should have the same number of templates. The first DNA template should give rise to the first peptide template and so on.
LibraryDesign objects are defined in the config file. 
## Data
During analysis, data is stored as a Data object instance. Data is just a container for individual samples, which are stored as SequencingSample objects. Any number of DNA sequences can be a sample in principle, but in practice, most of the time one sample = a single .fastq file.
SequencingSample objects have four public attributes: 

	SequencingSample.name: sample name (as a str)
    SequencingSample.D: a list of DNA sequences (can be set as None)
	SequencingSample.Q: a list of Q score sequences (can be set as None)
	SequencingSample.P: a list of peptide sequences (can be set as None)
	
These lists are stored as numpy arrays: 1D prior to calling FastqParser.transform() or FastqParser.translate(), and 2D arrays the entire time after that; shape: (number of entries, sequence length). Because the sequences for different reads may have a different length, arrays are padded to the longest sequence.\
\
The number of entries in each array is kept equal throughout the process unless one or more of the attributes are set to None. Although any given filtration routine [for example, FastqParser.q_score_filt()] acts on a single array [SequencingSample.Q in this example], the entries for all three arrays are discarded/kept as a result.\
\
Depending on how many and what kinds of templates are specified in LibraryDesign, any given entry in SequencingSample may in principle be compatible with several templates simultaneously. Figuring out what entry should be assigned to what kind of template is one of the primary objectives of the parser. Initially, [i.e. right after calling FastqParser.translate()] the parser deems every sequence to be compatible with every specified template. As filtration goes on, op by op, this compatibility is refined. I call the state of the assignment for a particular SequencingSample _sample’s internal state_. Some ops [for instance, FastqParser.fetch_at()] need to know exactly which template should be associated with which entry; if they find an entry that is compatible with multiple possible templates, they will “collapse” sample’s internal state, by choosing one compatible template and assigning everything else as incompatible. Refer to the list of ops below for details on which ops can collapse sample’s internal state. In general, these should be called after filtration ops.
## Ops

### FastqParser.fetch_fastq_from_dir()
    
    Fetch all .fastq files from the sequencing_data directory (as specified in config.py).
    Should be called as the first op in the workflow.
    
        Parameters:
                None
    
        Returns:
                Fetched Fastq data as an instance of Data

### FastqParser.fetch_gz_from_dir()
    
    Fetch all .fastq.gz files from the sequencing_data directory (as specified in config.py)
    Should be called as the first op in the workflow.
    
        Parameters:
                None
    
        Returns:
                Fetched Fastq data as an instance of Data

### FastqParser.revcom()
    
    For each sample in Data, get reverse complement of DNA sequences and 
    reverse sequences of the corresponding Q score. If used, should enqueued 
    right after the fetching op, and before any downstream ops.
    
        Parameters:
                None
    
        Returns:
                Transformed Data object holding reverse-complemented DNA and reversed Q score information

### FastqParser.transform()
    
    Deprecated in favor of using FastqParser.translate(). If used, should be called after 
    fetching the data and (optionally) running the FastqParser.revcom() op. Transforms 
    the data to a representation suitable for downstream ops.
    
        Parameters:
                None
    
        Returns:
                Transformed Data object

### FastqParser.translate(force_at_frame=None)
    
	For each sample in Data, perform in silico translation for DNA sequencing data. 
	The op will return data containing translated peptide lists. If used, should be 
    called after fetching the data and (optionally) running the FastqParser.revcom() op.
    
        Parameters:
                force_at_frame: if None, a regular ORF search will be performed. Regular ORF
                                search entails looking for a Shine-Dalgarno sequence upstream 
                                of an ATG codon (the exact 5’-UTR sequence signalling an 
                                ORF should be specified in config.py).
                                								
                                if not None, can take values of 0, 1 or 2. This will force-start
                                the translation at the specified frame regardless of the 
                                presence or absence of the SD sequence.
                                
                                For example:
                                DNA: TACGACTCACTATAGGGTTAACTTTAAGAAGGA
                   force_at_frame=0  ----------> 
                    force_at_frame=1  ---------->
                     force_at_frame=2  ---------->
					 
        Returns:
                Data object containing peptide sequence information

### FastqParser.len_filter(where=None, len_range=None)
    
    For each sample in Data, filter out sequences longer/shorter than the specified 
    library designs. Alternatively, a length range of sequences to take can be optionally 
    specified to filter out the entries (NGS reads) outside of this range.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
               len_range: either None (filtration will be done according to
                          the library design rules), or a list of two ints 
                          that specifies the length range to fetch.						  
					 
        Returns:
                Transformed Data object containing length-filtered data
				
### FastqParser.cr_filter(where=None, loc=None, tol=1)
    
    For each sample in Data, filter out sequences not containing intact constant
    regions. Entries (NGS reads) bearing constant regions with amino acids outside
	of the library design specification will be discarded.    
	
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify which constant regions 
                          the op should process. 

                     tol: int; specifies the maximum allowed number of mutations
                          constant region fetched with where/loc before the 
                          entry (NGS read) is discarded. For the library from above
                          
                seq:      ACDEF11133211AWVFRTQ12345YTPPK
             region:      [-0-][---1--][--2--][-3-][-4-]
        is_variable:      False  True   False True False
                          
                          calling cr_filter(where='pep', loc=[2], tol=1), will
                          discard all sequences containing more than 1 mutation
                          in the 'AWVFRTQ' region. Note that the insertions/deletions
                          in the constant region are not validated by the parser.					  
					 
        Returns:
                Transformed Data object containg entries with intact 
                constant regions
 				
### FastqParser.vr_filter(where=None, loc=None, sets=None)
    
    For each sample in Data, filter out sequences not containing intact variable 
    regions. Entries (NGS reads) bearing variable regions with amino acids outside
	of the library design specification will be discarded.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify which variable regions 
                          the op should process. 

                    sets: a list of ints; a list of monomer subsets to
                          check. For the library from above
                          
                seq:      ACDEF11133211AWVFRTQ12345YTPPK
             region:      [-0-][---1--][--2--][-3-][-4-]
        is_variable:      False  True   False True False
                          
                          there are five distinct variable amino acids:
                          1, 2, 3, 4, 5. The config file specifies which specific
                          amino acids are allowed for each of these numbers.
                          <vr_filter> op will make sure that each variable position
                          contains only the "allowed" monomers.					

                          vr_filter(where='pep', loc=[1], sets=[1, 3]) will make
                          sure that in region loc=1, variable amino acids 1 and 3
                          match the specification; variable amino acid 2 will not
                          be checked against in this example. Passing loc=[2] to
                          <vr_filter> op will raise an error, because it isn't a
                          variable region.
					 
        Returns:
                Transformed Data object containg entries with intact 
                variable regions
				
				
### FastqParser.filt_ambiguous(where=None)
    
    For each sample in Data, filter out sequences not containing intact ambiguous 
    tokens. For DNA, these are "N" nucleotides, which Illumina NGS routines occasionally
    assign during base calling. For peptides, these are any sequences containing
    amino acids outside of the translation table specification.	
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
        Returns:
                Transformed Data object containg entries without ambiguous
                tokens
				
### FastqParser.drop_data(where=None)
    
    For each sample in Data, delete datasets specified in 'where'. See documentation 
    on Data objects above for more information.
    
        Parameters:
                   where: 'dna', 'pep' or 'q' to specify which datasets 
                          should be dropped. 				
						  
        Returns:
                Transformed Data object without dropped datasets
				
### FastqParser.q_score_filt(minQ=None, loc=None)
    
    For each sample in Data, filter out sequences associated with Q scores below 
    the specified threshold minQ.
    
        Parameters:
                     loc: a list of ints to specify which regions 
                          the op should process. 

                    minQ: every Q score in the regions specified 
                          by loc should be greater or equal than 
						  this value; everything else will be discarded
                        						  
						  
        Returns:
                Transformed Data object
		
### FastqParser.fetch_at(where=None, loc=None)
    
    For each sample in Data, for a dataset specified by 'where', fetch the regions
    specified by 'loc'. Other regions are discarded. 
    
    Collapses sample's internal state (see above explanation for Data objects)
    
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify regions to be fetched 
						  
        Returns:
                Transformed Data object				

### FastqParser.unpad()
    
    For each sample in Data, unpads the D, Q, P arrays. For each array, removes 
	the columns where every value is a padding token. See documentation on Data 
    objects above for more information.

        Parameters:
                None	
						  
        Returns:
                Transformed Data object		
				
### FastqParser.len_summary(where=None, save_txt=False)
    
    For each sample in Data, compute the distribution of peptide/DNA sequence lengths
    (specified by 'where') and plot the resulting histogram in the parser output folder
    as specified by config.py. Optionally, the data can also be written to a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          						  
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots				
						  
        Returns:
                Data object (no transformation)
				
### FastqParser.convergence_summary(where=None)
    
    For each sample in Data, perform basic library convergence analysis on a sequence 
    level. Computes normalized Shannon entropy, and postition-wise sequence conservation. 
    Plots the results in the parser output folder as specified by config.py.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          						  							  
        Returns:
                Data object (no transformation)
				
				
### FastqParser.freq_summary(where=None, loc=None, save_txt=False)
    
    Perform basic library convergence analysis on a token level. For each sample in Data, 
    computes the frequency of each token in the dataset. Plots the results in the parser 
    output folder as specified by config.py. Optionally, the data can also be written to
    a txt file.
    
    Collapses sample's internal state (see above explanation for Data objects)
    	
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify regions to be analyzed

                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)

### FastqParser.q_summary(loc=None, save_txt=False)
    
    For each sample in Data, compute basic statistics of Q scores. 
	For each position in regions specified by 'loc', computes the mean and standard deviation
    of Q scores. Plots the results in the parser output folder as specified by config.py.
    Optionally, the data can also be written to a txt file.
    	
    Collapses sample's internal state (see above explanation for Data objects)
    	
        Parameters:					  
                     loc: a list of ints to specify regions to be analyzed

                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)	
				
### FastqParser.count_summary(where=None, top_n=None, fmt=None)
    
    For each sample in Data, counts the number of times each unique sequence is found in the
    dataset specified by 'where'. The results are written to a file in the parser output  
    folder as specified by config.py.
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                   top_n: if None, full summary will be created. If
                          an int is passed, only top_n sequences (by count)
                          will be written to a file.

                     fmt: the format of the output file. Supported values are
                          'csv' and 'fasta'.					 
                          						  							  
        Returns:
                Data object (no transformation)
				
### FastqParser.save(where=None, fmt=None)
    
    For each sample in Data, save the dataset specified by 'where'. The results are written 
    to a file in the parser output folder as specified by config.py.
    
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                     fmt: the format of the output file. Supported values are
                          'npy', 'fasta' and 'csv'					 
                          						  							  
        Returns:
                Data object (no transformation)

# Known issues
1. After FastqParser.fetch_at() is called, the assignment between fetched sequence regions and LibaryDesign gets broken. It means that no other ops that takes 'loc' as a keyword can be called after FastqParser.fetch_at() [ops that don't take 'loc' as a keyword are ok]. To be fixed.
2. FastqParser.q_summary() collapses sample's internal state, which can often be inconvenient. To be fixed.