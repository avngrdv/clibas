# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 00:19:28 2022
@author: a_vin
"""

import config
from clibas.dispatchers import Dispatcher
from clibas.pipelines import Pipeline
from clibas.parsers import FastqParser
#from clibas.baseclasses import Handler
from clibas.dataanalysis import HDBUMAP, DataAnalysisTools

import numpy as np

# handlers = (Pipeline, FastqParser, HDBUMAP)
# pip, par, hdu = Dispatcher(config).dispatch_handlers(handlers)

handlers = (Pipeline, FastqParser, DataAnalysisTools)
pip, par, dat = Dispatcher(config).dispatch_handlers(handlers)

# pip.enque([
#             par.fetch_gz_from_dir(), 
#             par.translate(stop_readthrough=False),         
#             par.len_filter(where='pep'),      
#             par.cr_filter(where='pep', loc=[1], tol=3),
#             par.vr_filter(where='pep', loc=[0], sets=[1, 2, 3]),
#             par.filt_ambiguous(where='pep'),
#             par.q_score_filt(minQ=30, loc=[1]),
#             par.fetch_at(where='pep', loc=[0]),
#             par.unpad(),
#             par.save(where='pep', fmt='npy')
#           ])

# data = pip.run()

F = np.load('DENSE_Morgan_F_r=4_LazDEF.npy')    
# # #F = np.load('DENSE_Morgan_F_r=4_LazBF.npy')

# #bCD
# # arr1 = np.load('Yunxiang1_S1_L001_R1_001.fastq_pep.npy').astype('<U1')
# # arr2 = np.load('Yunxiang2_S2_L001_R1_001.fastq_pep.npy').astype('<U1')
# # arr3 = np.load('Yunxiang3_S3_L001_R1_001.fastq_pep.npy').astype('<U1')
# # arr4 = np.load('Yunxiang4_S4_L001_R1_001.fastq_pep.npy').astype('<U1')
# # arr5 = np.load('Yunxiang5_S5_L001_R1_001.fastq_pep.npy').astype('<U1')
# # arr6 = np.load('Yunxiang6_S6_L001_R1_001.fastq_pep.npy').astype('<U1')

# # #Harry
# arr1 = np.load('Pf4_S46_L001_R1_001.fastq_pep.npy').astype('<U1')
# arr2 = np.load('Pf5_S47_L001_R1_001.fastq_pep.npy').astype('<U1')
# arr3 = np.load('Pf6_S48_L001_R1_001.fastq_pep.npy').astype('<U1')
# arr4 = np.load('Pf7_S49_L001_R1_001.fastq_pep.npy').astype('<U1')

# # #TNIK
# arr1 = np.load('tnik_r3_pre_for_p_matrix.npy').astype('<U1')
# arr2 = np.load('tnik_r4_pre_for_p_matrix.npy').astype('<U1')
# arr3 = np.load('tnik_r5_pre_for_p_matrix.npy').astype('<U1')

# # # #Heinis
# # arr1 = np.load('FXIa_Round2_Linker7_pep.npy').astype('<U1')
# # arr2 = np.load('FXIa_Round3_Linker7_pep.npy').astype('<U1')
# # arr3 = np.load('FXIa_Round5_Linker7_0pSIF_pep.npy').astype('<U1')
# # arr4 = np.load('FXIa_Round5_Linker7_1pSIF_pep.npy').astype('<U1')
# # arr5 = np.load('FXIa_Round5_Linker7_10pSIF_pep.npy').astype('<U1')
# # arr6 = np.load('FXIa_Round6_Linker7_0pSIF_pep.npy').astype('<U1')
# # arr7 = np.load('FXIa_Round6_Linker7_1pSIF_pep.npy').astype('<U1')
# # arr8 = np.load('FXIa_Round6_Linker7_10pSIF_pep.npy').astype('<U1')

# # # # arr[arr == 'Z'] = 'C'
# # # # arr = np.load('Pf5_S47_L001_R1_001.fastq_pep.npy').astype('<U1')

import pandas as pd
df = pd.read_csv('peptides.csv')
arr = list(df['Peptide'])

# from clibas.misc import sample_from_template
from clibas.datatypes import Data, AnalysisSample


# # titration = [0.005, 0.01, 0.02, 0.05, 0.1, 0.3]
# # tot = 3000
# # samples = []

# # for x in titration:
# #     arr1 = sample_from_template('XXXXXXXXXXXXXX', int(tot * (1 - 2 * x)), hdu.constants.aas)
# #     arr2 = sample_from_template('VIXXARXXGIXLXX', int(tot * x), hdu.constants.aas)
# #     arr3 = sample_from_template('XXNEXXXVEXXRXY', int(tot * x), hdu.constants.aas)
# #     arr = np.vstack((arr1, arr2, arr3))
# #     samples.append(AnalysisSample(X=arr, name=f'motif fraction {x}')) 
                    
# # data = Data(samples)

# data = Data([
#               AnalysisSample(X=arr1, name='tnik r3'),
#               AnalysisSample(X=arr2, name='tnik r4'),
#               AnalysisSample(X=arr3, name='tnik r5'),
#               # AnalysisSample(X=arr4, name='harry r7'),
#               # AnalysisSample(X=arr5, name='FXIa phage round 5 10SIF'),
#               # AnalysisSample(X=arr6, name='FXIa phage round 6 0SIF'),
#               # AnalysisSample(X=arr7, name='FXIa phage round 6 1SIF'),
#               # AnalysisSample(X=arr8, name='FXIa phage round 6 10SIF'),             
#             ])

# hdu.analyze(top_n=3000,
#             where='pep',
#             F=F,
#             plot_freqs=True,
#             cluster_fasta=False,
#             single_manifold=False)(data)



# # (X, Y, C, labels, cluster_scores) = hdu.embed_and_cluster(data[-1].X, 
# #                                                           top_n=3000,
# #                                                           alphabet=hdu.constants.aas, 
# #                                                           F=F
# #                                                          )

# # import clibas.plotters as P
# # from clibas.misc import get_freqs
# # # P.Analysis.UMAP_HDBSCAN(Y, C, labels)


# # # cluster = 3
# # aas = hdu.constants.aas
# # freq = get_freqs(data[3].X, aas)
# # P.SequencingData.tokenwise_frequency(freq, aas)

# # arr = hdu._recast_data(arr)
# # freq = get_freqs(arr[0].X, aas)
# # P.SequencingData.tokenwise_frequency(freq, aas)










