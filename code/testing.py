# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 00:19:28 2022
@author: a_vin
"""

import config
from clibas.dispatchers import Dispatcher
from clibas.pipelines import Pipeline
from clibas.parsers import FastqParser
from clibas.dataanalysis import DataAnalysisTools

handlers = (Pipeline, FastqParser, DataAnalysisTools)
pip, par, dat = Dispatcher(config).dispatch_handlers(handlers)

pip.enque([
            par.fetch_gz_from_dir(), 
            par.translate(stop_readthrough=False),         
            dat.length_analysis(where='dna', save_txt=True),
            dat.length_analysis(where='pep', save_txt=True),
            par.len_filter(where='pep'),
            dat.sequence_convergence_analysis(where='dna'),
            dat.sequence_convergence_analysis(where='pep'),            
            par.cr_filter(where='pep', loc=[1], tol=3),
            par.vr_filter(where='pep', loc=[0], sets=[1, 2, 3]),
            par.filt_ambiguous(where='pep'),
            par.q_score_filt(minQ=20, loc=[1]),
            dat.token_convergence_analysis(where='pep', loc=None, save_txt=True),
            par.fetch_at(where='pep', loc=[0]),
            dat.token_convergence_analysis(where='pep', loc=[0], save_txt=True),
            dat.token_convergence_analysis(where='dna', loc=[1], save_txt=True),
            par.fastq_count_summary(where='pep', top_n=500, fmt='csv'),
            par.fastq_count_summary(where='pep', top_n=500, fmt='fasta'),
            par.fastq_count_summary(where='dna', top_n=500, fmt='csv'),
            par.fastq_count_summary(where='dna', top_n=500, fmt='fasta'),      
            dat.q_score_analysis(loc=None, save_txt=True),
            par.fetch_at(where='dna', loc=[0]),
            dat.q_score_analysis(loc=[0], save_txt=True),
            par.library_design_match_analysis(where='pep'),
            par.unpad(),
            par.save(where='pep', fmt='npy'),
            dat.umap_hdbscan_analysis(where='pep',
                                      top_n=3000,
                                      F='pep_ECFP4', 
                                      cluster_fasta=True,
                                      single_manifold=True,
                                      return_modified=True)
          ])

data = pip.run()