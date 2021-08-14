# -*- coding: utf-8 -*-
"""
Created on Sat Aug 8 20:28:10 2021
@author: Alex Vinogradov
"""
if __name__ == '__main__':
    
    #import prerequisities
    from utils.ProcessHandlers import Pipeline, FastqParser
    from utils.Dispatcher import Dispatcher
    
    #config file holds the information about library designs and
    #other parser instructions (where to look for data, where to save results etc)
    import config
    
    #initialize a dispatcher object; dispatcher is strictly speaking
    #not necessary, but it simplifies initialization of data handlers
    dispatcher = Dispatcher(config)    
    
    #a list of handlers to initialize; pipeline should always be included
    #if NGS data parsing is the goal, FastqParser will do most of the work
    handlers = (Pipeline, FastqParser)
    
    #initialize the handlers
    pip, par = dispatcher.dispatch_handlers(handlers)
    
    '''
    FastqParser op keywords: 
        
    where - Which dataset should be used to perform filtration/analysis. 
            Either 'pep' or 'dna' for most operations.

    loc -  A list of pointers to library design regions.
           For example, for the following library sequence:
               
                    seq:      ACDEF11133211AWVFRTQ12345YTPPK
                    loc:      [-0-][---1--][--2--][-3-][-4-]
            is_variable:      False  True   False True False
    
           To perform a constant region check, (i.e. to make sure
           that every sequence in the sample contains an intact 
           constant region subsequence as specified by the library design),
           loc can be specified as [0, 2, 4].  Or, if we only care about
           the 'AWVFRTQ' subsequence being intact, loc can be specified as [2].
          
    tol -  Integer; maximum tolerated number of mutations for
           a constant region. Going with the example from above,
           typing cr_filter(where='pep', loc=[2], tol=1), will
           discard all sequences containing more than 1 mutation
           in the 'AWVFRTQ' region. Note that the insertions/deletions
           in the constant region are not tolerated by the parser.
          
    sets - A list of monomer subsets to check. For the example
           above, there are five distinct variable amino acids:
           1, 2, 3, 4, 5. The config file specifies which *specific* 
           amino acids are allowed for each of these numbers.
           <vr_filter> op will make sure that each variable position
           contains only the "allowed" monomers.
           
           vr_filter(where='pep', loc=[1], sets=[1, 3]) will make
           sure that in region loc=1, variable amino acids 1 and 3
           match the specification; variable amino acid 2 will not
           be checked in this example. Passing loc=[2] to <vr_filter> 
           op will raise an error, because it isn't a variable region.
           
    minQ - Minimal Q score value allowed for a DNA sequence
           in the region at loc. A parameter for <q_score_filt> 
           
    save_txt  -  False/True; for ops performing data analysis, i.e.
                 <len_summary>, <freq_summary>, <q_summary>
                 
    fmt       -  A param for ops performing data analysis, namely:
                 <save>: fmt is any of 'fasta', 'csv', 'npy'
                 <count_summary>: fmt is any of 'fasta', 'csv'
 
    force_at_frame - an optional parameter for the <translate> op.
                     If a DNA read begins in the middle of an ORF
                     and thus is missing the ORF promoter region
                     (i.e. SD sequence + at ATG codon downstream),
                     translation function will return nothing.
                     
                     In such cases, force_at_frame can be set to
                     any of (0, 1, 2) to force-start the translation
                     at the specified frame.
                     
                     dna: CGACTCACTATAGGGTTAACTTTAAGAAGGAGAT
           force_at_frame=0--------------->
            force_at_frame=1-------------->   
             force_at_frame=2------------->
    '''

    
    #enqueue the list of ops to run
    #note that at this stage no data processing will take place,
    #but the validity of passed commands will be validated.
    pip.enque([
                par.fetch_gz_from_dir(), 
                par.translate(),         
                par.len_summary(where='dna', save_txt=True),
                par.len_summary(where='pep', save_txt=True),
                par.len_filter(where='pep'),          
                par.convergence_summary(where='dna'),
                par.convergence_summary(where='pep'),            
                par.cr_filter(where='pep', loc=[1], tol=2),
                par.vr_filter(where='pep', loc=[0], sets=[1, 2, 3]),
                par.q_summary(loc=[0, 1], save_txt=True),
                par.filt_ambiguous(where='pep'),
                par.q_score_filt(minQ=20, loc=[0]),
                par.freq_summary(where='pep', loc=[0], save_txt=True),
                par.freq_summary(where='dna', loc=[0, 1], save_txt=True),
                par.fetch_at(where='pep', loc=[0]),
                par.count_summary(where='pep', top_n=50, fmt='csv'),
                par.unpad(),
                par.save(where='pep', fmt='fasta'),
             ])
    
        
    #this will execute the pipeline
    #if save_summary=True, summary will be saved in the ../logs folder
    #as specified in the config file
    data = pip.run(save_summary=True)
        
        
    
    
    
    
    

















