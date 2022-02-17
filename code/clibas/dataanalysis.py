# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 22:42:26 2022
@author: Alex Vinogradov
"""

import os, inspect
import numpy as np
import pandas as pd

import clibas.plotters as Plotter
from clibas.baseclasses import Handler
from clibas.datatypes import Data, AnalysisSample

class DataAnalysisTools(Handler):
    '''
    Crafted wrappers around various data analysis tools to simplify their
    calling during the pipeline creation process.
    '''
    
    def __init__(self, *args):
        super(DataAnalysisTools, self).__init__(*args)
    
        self._validate_designs()
        self._validate_constants()
        self._on_completion()
        return

    def __repr__(self):
        return '<DataAnalysisTools object>'
    
    def length_analysis(self, where=None, save_txt=False):
        
        self._where_check(where)
        la = LengthAnalysis(self.__dict__)
        op = la.len_summary(self, where=where, save_txt=save_txt)
        return op
    
    def sequence_convergence_analysis(self, where=None, alphabet=None):

        self._where_check(where)
        alphabet = self._infer_alphabet(where, alphabet)
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.sequence_level_convergence(self, where=where, alphabet=alphabet)
        return op
        
    def token_convergence_analysis(self, where=None, 
                                   loc=None, 
                                   alphabet=None,
                                   save_txt=False):
    
        self._where_check(where)
        design = self._infer_design(where)
        alphabet = self._infer_alphabet(where, alphabet)
        if loc is not None:
            self._loc_check(loc, design)
        
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.token_level_convergence(self, 
                                        where=where, 
                                        loc=loc, 
                                        alphabet=alphabet,
                                        save_txt=save_txt
                                       )
        return op
        
    def umap_hdbscan_analysis(self, 
                              top_n=None,
                              where=None,
                              F=None, 
                              plot_freqs=False,
                              cluster_fasta=False,
                              alphabet=None,
                              single_manifold=False):

        hdu = HDBUMAP(self.__dict__)
        op = hdu.analyze(top_n=top_n, where=where,
                         F=F, plot_freqs=plot_freqs,
                         cluster_fasta=cluster_fasta,
                         alphabet=alphabet, single_manifold=single_manifold
                        )
        
        return op
        

    def q_summary(self, loc=None, save_txt=False):
        '''
        For each sample in Data, compute some basic Q score statistics.
    	For each position in regions specified by 'loc', computes the mean and standard deviation
        of Q scores. Plots the results in the parser output folder as specified by config.py.
        Optionally, the data can also be written to a txt file.
        	    	
        Parameters:					  
                     loc: a list of ints to specify regions to be analyzed;
                          in this case, the op will collapse sample's internal
                          state (see explanation for Data objects)
                          
                          OR
                          
                          'all': to get the same statistics over the entire Q 
                                 score arrays. in this case, the op will NOT 
                                 collapse sample's internal state
                          
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)	
        '''
        
        if loc != 'all':
            self._loc_check(loc, self.D_design)
            
        def q_score_summary(data):
            
            #here recasting really doesn't make sense
            self._prepare_destinations(data, self.dirs.analysis_out)
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                
                if loc == 'all':
                    relevant_arr = sample.Q.astype(np.float32)
                    
                else:
                    arr = sample.Q
                    if not sample._is_collapsed:
                        msg = f"<q_score_summary> routine will collapse sample {sample.name}'s internal state"
                        self.logger.info(msg)
                        sample._collapse_internal_state()
                                    
                    maxlen = self._find_max_len(self.D_design, loc)
                    #iterate over templates and append all of the relevant arr views to this array
                    #relevant view: masked (row/columnwise) arr
                    relevant_arr = []
                    
                    for i,template in enumerate(self.D_design):
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)                
                        
                        arr_view = np.zeros((np.sum(row_mask), maxlen), dtype=np.float32)
                        arr_view[:,:len(col_mask)] = arr[row_mask][:,col_mask]
                        relevant_arr.append(arr_view)
    
                    #assemble into a single array
                    relevant_arr = np.vstack(relevant_arr)
                    
                #mask out pads (0) as nans for nanmean/nanstd statistics
                relevant_arr[relevant_arr == 0] = np.nan
                
                #get the stats; plot
                q_mean = np.nanmean(relevant_arr, axis=0)
                q_std = np.nanstd(relevant_arr, axis=0)
    
                destination = os.path.join(self.dirs.parser_out, sample.name)
                if loc == 'all':
                    nloc = 'overall'
                    fname = f'{sample.name}_overall_q_score_summary'
                else:
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_reg{nloc}_q_score_summary'
                    
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.Q_score_summary(q_mean, q_std, nloc, basename)  
    
                if save_txt:
                    q = np.vstack((q_mean, q_std))
                    np.savetxt(basename + '.csv',
                               q.T,
                               delimiter=',',
                               header='Q mean, Q std')  
                
            return data
        return q_score_summary
    
    def count_summary(self, where=None, top_n=None, fmt=None):
        '''
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
        '''
        self._where_check(where)
        
        if fmt not in ('csv', 'fasta'):
            msg = f"<count_summary> routine received invalid fmt argument. Acceted any of ('csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)
            
        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f'<count_summary> routine expected to receive parameter top_n as as int; received: {type(top_n)}'
                self.logger.error(msg)  
                raise ValueError(msg)
    
        def _writer(sample, og_ind, counts, fmt, path):
            if fmt == 'csv':           
                df = pd.DataFrame(columns=['Peptide', f'{where} count', 'DNA'])
                df['Peptide'] = [''.join(x) for x in sample.P[og_ind]]
                df['DNA'] = [''.join(x) for x in sample.D[og_ind]]
                df[f'{where} count'] = counts
                df.to_csv(path + '.csv', sep=',')
                
            if fmt == 'fasta':
                
                arr = sample[where][og_ind]
                arr_1d = [''.join(x) for x in arr]    
                
                with open(path + '.fasta', 'w') as f:
                    for i,seq in enumerate(arr_1d):
                        f.write(f'>seq_{i+1}_count_{counts[i]}\n')  
                        f.write(f'{seq}\n')
                return 
        
        def full_count_summary(data):
            self._prepare_destinations(data, self.dirs.analysis_out)
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]            
                
                #count entries in the array
                unique, og_ind, counts = np.unique(arr, axis=0, 
                                                   return_counts=True,
                                                   return_index=True)                
                
                #if top_n is unset, ind array will index every entry in the sample
                ind = np.argsort(counts)[::-1][:top_n]
    
                og_ind = og_ind[ind]
                counts = counts[ind]
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_count_summary'
                path = os.path.join(destination, fname)
                
                _writer(sample, og_ind, counts, fmt, path)
                                
            return data
        return full_count_summary
    
    def template_summary(self, where=None):
        '''
        For each sample in Data, compute the number of matches between the dataset 
        specified by 'where' and the corresponding library templates. The results 
        are written to a file in the parser output folder as specified by config.py.
        
        In other words, summarize where dataset sequences come from (from which
        libraries). The op could also be called "_internal_state_summary"
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on                        						  							  
        Returns:
                Data object (no transformation)
        '''        
        
        self._where_check(where)
        design = self._infer_design(where)
            
        def template_breakdown(data):       
            #summarize straight into a pandas dataframe
            sample_names = [sample.name for sample in data]
            templates = [template.lib_seq for template in design]
        
            #all this op is: axis=0-wide sum of the internal states
            import pandas as pd
            df = pd.DataFrame(index=sample_names, columns=templates)
            
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                df.loc[sample.name] = np.sum(sample._internal_state, axis=0)
                
            fname = f'{self.exp_name}_by_template_breakdown.csv'
            path = os.path.join(self.dirs.logs, fname)            
            df.to_csv(path + '.csv', sep=',')
    
            return data
        return template_breakdown
  
    
class HDBUMAP(Handler):
    
    def __init__(self, *args):
        super(HDBUMAP, self).__init__(*args)
        self._validate_imports()
        return

    def _validate_imports(self):
        try:
            import umap
            import hdbscan
            
            self.umap = umap
            self.hdbscan = hdbscan
        except:
            msg = 'Failed to import the libraries necessary for HDBSCAN/UMAP analysis. . .'
            self.logger.error(msg)
            raise ImportError(msg)
        return

    def __repr__(self):
        return '<HDBSCAN UMAP analysis object>'    
    
    def _preprocess_arr(self, arr, top_n=None, alphabet=None, F=None):
        
        from clibas.datapreprocessors import DataPreprocessor
        from clibas.misc import sorted_count, arr_purity
        
        #this casting to AnalysisSample is done to ensure that a 2D array
        #in a proper format is passed downstream
        arr = AnalysisSample(X=arr)        
        X, C = sorted_count(arr.X, top_n=top_n)
        p_arr = self._recast_data(X)
        
        #featurize the array; everything is done in memory, so the array has to fit
        pre = DataPreprocessor(dict(F=F))
        p_arr = pre.int_repr(alphabet=alphabet)(p_arr)
        p_arr = pre.featurize_X(alphabet=alphabet,
                                reshape=False,
                                repad=False)(p_arr)

        #setup umap/clustering hyperparams right away to save computation
        purity = arr_purity(arr, alphabet)
        arr_size = X.shape[0]
        self._infer_hyperparameters(arr_size, purity)

        return X, C, p_arr[0]['X']
    
    def _initialize_umap(self):
        
        if not hasattr(self, 'n_neighbors'):
            self.n_neighbors = 15
        
        if not hasattr(self, 'min_dist'):
            self.min_dist = 0.1
            
        self.embedder = self.umap.UMAP(n_neighbors=self.n_neighbors,
                                       min_dist=self.min_dist,
                                       n_components=2,
                                       spread=1.0
                                      )
        return
      
    def _infer_hyperparameters(self, size, purity):
        '''
        use heuristic parametric curves to infer hyperparameters
        parameter curves are a function of dataset purity: datasets with
        low purity need to be repsented on a more "local" manifold is the
        gist.
        
        arr: array to be embedded'''
        
        def hbscan_param_curve(heuristic, size, purity):
            
            return np.ceil(np.divide(heuristic * np.log2(size), 
                                     -np.log10(purity))
                          )
          
        #4PL params are heuristic
        from clibas.misc import logistic_4_param
        self.min_dist = logistic_4_param(purity, 0.01, 0.65, 20, 0.25)
                
        if size < 200:
            self.n_neighbors = 5
        
        else:
            n_neighbors_min = max(0.001 * size, 5)
            n_neighbors_max = min(0.04 * size, 95)
            self.n_neighbors = np.ceil(logistic_4_param(purity,
                                                        n_neighbors_min,
                                                        n_neighbors_max,
                                                        15, 0.35
                                                       )
                                       ).astype(np.int)

        print(f'purity: {purity:.3f}, min_dist: {self.min_dist:.3f}, n_neighbors: {self.n_neighbors}')

        #for numerical stability
        purity += 0.001
           
        #HBDSCAN params also need to be approximated using array purity as a proxy
        min_cluster_heur = np.logspace(-4, 4, num=8, base=1.4)
        min_samples_heur= np.logspace(-10, 0, num=8, base=1.4)
        
        min_cluster_size = hbscan_param_curve(min_cluster_heur, size, purity)
        min_cluster_size[min_cluster_size < 2] = 2
        self.min_cluster_size = np.unique(min_cluster_size).astype(np.int)
        
        min_samples = hbscan_param_curve(min_samples_heur, size, purity)
        min_samples[min_samples < 1] = 1
        self.min_samples = np.unique(min_samples).astype(np.int)
        
        print(f'min_cluster_size: {self.min_cluster_size}')
        print(f'min_samples: {self.min_samples}')
        return
    
    def _transform_scale(self, pX):
        from sklearn.preprocessing import StandardScaler        
        Y = self.embedder.transform(pX)
        return StandardScaler().fit_transform(Y)

    def _cluster(self, X, Y, alphabet=None):
        
        #Y: scaled umap embedding
        if not hasattr(self, 'min_cluster_size'):
            self.min_cluster_size = np.array([5])
            #as array to make it iterable
        
        if not hasattr(self, 'min_samples'):
            self.min_samples = np.array([3])
            
        from clibas.misc import naive_clustering_score
        cluster_scores = np.zeros((self.min_cluster_size.size,
                                   self.min_samples.size
                                 ))
    
        #scout for the best clustering hyperparameters
        for i,c in enumerate(self.min_cluster_size):
            for j,s in enumerate(self.min_samples):
                
                hdb = self.hdbscan.HDBSCAN(
                                           cluster_selection_epsilon=0.01,
                                           min_cluster_size=int(c),
                                           min_samples=int(s),
                                           cluster_selection_method = 'eom'
                                          )
                
                hdb.fit(Y)
                cluster_scores[i, j] = naive_clustering_score(X,
                                                              hdb.labels_,
                                                              alphabet=alphabet,
                                                              return_mean=True
                                                             ) 
        
        #use the best hyperparamaters to return the actual labels
        i, j = np.where(cluster_scores == cluster_scores.max())
        hdb = self.hdbscan.HDBSCAN(
                                   cluster_selection_epsilon=0.01,
                                   min_cluster_size=int(self.min_cluster_size[i[0]]),
                                   min_samples=int(self.min_samples[j[0]]),
                                   cluster_selection_method = 'eom'
                                  )

        hdb.fit(Y)
        return hdb.labels_, cluster_scores
    
    def embed_and_cluster(self, arr, top_n=None, alphabet=None, F=None):
        #TODO: is there a modular way to run all of the hyperparameter checks
        #like in the analyze op?
        
        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        self._initialize_umap()
        self.embedder.fit(pX)        
        Y = self._transform_scale(pX)
        labels, cluster_scores = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_scores)

    def transform_and_cluster(self, arr, top_n=None, alphabet=None, F=None):

        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        Y = self._transform_scale(pX)
        labels, cluster_scores = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_scores)

    def cluster_summary(self, X, labels, alphabet=None):
        
        from clibas.misc import naive_clustering_score, arr_purity
        cluster_scores = naive_clustering_score(X, 
                                                labels,
                                                alphabet=alphabet, 
                                                return_mean=False
                                               )

        cluster_sizes = [X[labels == i].shape[0] 
                         for i in np.unique(labels)
                        ]
        
        #inefficient computation, but i can't imagine this ever being the bottleneck
        cluster_purities = [arr_purity(X[labels ==i], alphabet=alphabet) 
                            for i in np.unique(labels)
                           ]
    
        d = {'Cluster purity': cluster_purities,
             'Cluster size': cluster_sizes,
             'Cluster number': np.unique(labels),
             'Cluster score': cluster_scores,
            }
        
        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(['Cluster score'], ascending=False)    
        return df

    def entry_count_cluster(self, X, C, labels):
    
        d = {'Entry': np.array([''.join(x) for x in X]),
             'Count': C,
             'Cluster': labels}         
        
        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(['Cluster', 'Count'], ascending=False)
        return df

    def analyze(self,
                top_n=None,
                where=None,
                F=None, 
                plot_freqs=False,
                cluster_fasta=False,
                alphabet=None,
                single_manifold=False):
                
        #TODO: document
        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f'<HDBUMAP.analyze> op expected to receive parameter top_n as int; received: {type(top_n)}'
                self.logger.error(msg)  
                raise ValueError(msg)
        else:
            msg = '<HDBUMAP.analyze>: top_n param is unspecified. Op may take a long time to finish if the datasets are large. . .'
            self.logger.warning(msg)
            
        if not isinstance(cluster_fasta, bool):
            msg = f'<HDBUMAP.analyze> op expected to receive parameter cluster_fasta as bool; received: {type(cluster_fasta)}'
            self.logger.error(msg)  
            raise ValueError(msg)
            
        if not isinstance(plot_freqs, bool):
            msg = f'<HDBUMAP.analyze> op expected to receive parameter plot_freqs as bool; received: {type(plot_freqs)}'
            self.logger.error(msg)  
            raise ValueError(msg)        
        
        if F is not None:
            if not isinstance(F, np.ndarray):
                msg = f'<HDBUMAP.analyze>: op expected to receive parameter F as numpy array; received: {type(F)}'
                self.logger.error(msg)
                raise ValueError(msg)

            else:
                if not F.ndim == 2:
                    msg = f'<HDBUMAP.analyze>: op expected to receive parameter F as a two-dimensional array; received: ndims={F.ndim}'
                    self.logger.error(msg)
                    raise ValueError(msg) 
                    
        if alphabet is not None:
            if not isinstance(alphabet, (tuple, list, np.ndarray)):
                msg = f'<HDBUMAP.analyze> op expected to receive parameter alphabet as iterable (list, tuple, ndarray); received: {type(alphabet)}'
                self.logger.error(msg)  
                raise ValueError(msg)            
                
        else:
            alphabet = self._infer_alphabet(where)
        
        if not alphabet:
            msg = '<HDBUMAP.analyze> could not infer token alphabet. . .'
            self.logger.error(msg)  
            raise ValueError(msg)  
        
        if not isinstance(single_manifold, bool):
            msg = f'<HDBUMAP.analyze> op expected to receive parameter single_manifold as bool; received: {type(single_manifold)}'
            self.logger.error(msg)  
            raise ValueError(msg)
                
        def _cluster_writer(df, destination, sname):
            
            clusters = np.unique(df['Cluster'])
            destination = os.path.join(destination, 'fasta_files')
            
            if not os.path.isdir(destination):
                os.makedirs(destination)               
            
            for cluster in clusters:
            
                fname = f'{sname}_{where}_cluster_{cluster}.fasta'
                path = os.path.join(destination, fname)
                
                with open(path, 'w') as f:
                    for entry in df[df['Cluster'] == cluster].iterrows():
                        f.write(f'>cluster_{cluster}_count_{entry[1]["Count"]}\n')  
                        f.write(f'{entry[1]["Entry"]}\n')
            return           
        
        def _dump_freqs(X, labels, alphabet, destination, sname):
            
            #TODO: i think this should just call to the freq_analysis class, innit
            from clibas.misc import get_freqs
            clusters = np.unique(labels)
            destination = os.path.join(destination, 'consensus_plots')
            
            if not os.path.isdir(destination):
                os.makedirs(destination)     
            
            for cluster in clusters:
                
                freq = get_freqs(X[labels == cluster], alphabet)
                basename = os.path.join(destination, f'{sname}_cluster_{cluster}_consensus')
                Plotter.SequencingData.tokenwise_frequency(freq, alphabet,
                                                           loc=None,
                                                           where=f'cluster={cluster}',
                                                           basename=basename)

            return
        
        def full_analysis(data):
            
            #make sure any input is ok -> recast everything to Data(dtype=AnalysisSample)
            recast_data = self._recast_data(data, where=where)
            func = self.embed_and_cluster
            
            if single_manifold:
                stacked, recast_data = self._restack_and_repad(recast_data)
                
                _, _, pX = self._preprocess_arr(stacked, 
                                                top_n=top_n,
                                                alphabet=alphabet,
                                                F=F
                                               )
                self._initialize_umap()
                self.embedder.fit(pX)
                func = self.transform_and_cluster
                
            for sample in recast_data:
            
                X, Y, C, labels, scores = func(sample.X,
                                               top_n=top_n, 
                                               alphabet=alphabet,
                                               F=F
                                              )                
                
                destination = os.path.join(self.dirs.analysis_out, 
                                           sample.name, 
                                           'UMAP_HDBSCAN_analysis')
                
                if not os.path.isdir(destination):
                    os.makedirs(destination)       

                #plot the results
                basename = os.path.join(destination, 
                                        f'{sample.name}_{where}_UMAP_HDBSCAN')
                
                Plotter.Analysis.UMAP_HDBSCAN(Y, labels,
                                              C=C,
                                              sample_name=sample.name,
                                              basename=basename,
                                              show_annotations=True
                                             )

                #save cluster summary
                df = self.cluster_summary(X, labels, alphabet=alphabet)
                fname = f'{sample.name}_{where}_clustering_summary.csv'
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=',', index=False)                
            
                #make + save an entry-count-label summary
                df = self.entry_count_cluster(X, C, labels)
                fname = f'{sample.name}_{where}_entry_count_cluster.csv'
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=',', index=False)                     
                
                #plot clustering hyperparams
                basename = os.path.join(destination, 
                                        f'{sample.name}_{where}_clustering_optimization')                
                
                Plotter.Analysis.ClusteringHyperParams(self.min_cluster_size,
                                                       self.min_samples,
                                                       scores,
                                                       sample_name=sample.name,
                                                       basename=basename
                                                       )                
                
                if cluster_fasta:
                    _cluster_writer(df, destination, sample.name)
                
                if plot_freqs:
                    _dump_freqs(X, labels, alphabet, destination, sample.name)
            
            return data
        return full_analysis
  


class LengthAnalysis(Handler):
    
    def __init__(self, *args):
        super(LengthAnalysis, self).__init__(*args)


    def len_summary(self, where=None, save_txt=False):
        '''
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
        '''
        
        def length_summary(data):
            
            recast_data = self._recast(data, where=where)
            self._prepare_destinations(data, self.dirs.analysis_out)        
            for sample in recast_data:
                
                arr = sample.X
                L = self._L_summary(arr)            
                L, counts = np.unique(L, return_counts=True)
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_L_distribution'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.L_distribution(L, counts, 
                                                      where=where, 
                                                      basename=basename
                                                     )
                
                if save_txt:
                    np.savetxt(basename + '.csv',
                               np.array((L, counts)).T,
                               delimiter=',', 
                               header='Seq length,Count')
            
            return data
        return length_summary    
    
class ConvergenceAnalys(Handler):
    
    def __init__(self, *args):
        super(LengthAnalysis, self).__init__(*args)
        
        from clibas.misc import (
                                 shannon_entropy,
                                 get_freqs,
                                 positional_conservation
                                )

        self.shannon_entropy = shannon_entropy
        self.get_freqs = get_freqs
        self.positional_conservation = positional_conservation
        return

    def sequence_level_convergence(self, where=None, alphabet=None):
        '''
        For each sample in Data, perform basic library convergence analysis on 
        a sequence level. Computes normalized Shannon entropy and postition-wise 
        sequence conservation. Plots the results in the parser output folder as 
        specified by config.py.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          
                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        
        def sequence_level_convergence_summary(data):
    
            recast_data = self._recast(data, where)
            self._prepare_destinations(data, self.dirs.analysis_out)
            
            for sample in recast_data:
    
                self._transform_check(sample, inspect.stack()[0][3])                
                arr = sample.X
                
                shannon, counts = self.shannon_entropy(arr, norm=True)
                freq = self.get_freqs(arr, alphabet)
                seq_conservation = self.positional_conservation(freq)
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_library_convergence'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.dataset_convergence(counts, shannon, where, basename)                
                
                fname = f'{sample.name}_{where}_sequence_conservation'
                basename = os.path.join(destination, fname)                
                Plotter.SequencingData.conservation(seq_conservation, where, basename)
                
            return data
        return sequence_level_convergence_summary


    def token_level_convergence(self, where=None, loc=None, alphabet=None, save_txt=False):
        '''
        Perform basic library convergence analysis at a token level. For each sample in Data, 
        computes the frequency of each token in the dataset. Plots the results in the parser 
        output folder as specified by config.py. Optionally, the data can also be written to
        a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
    						  
                     loc: if not None: a list of ints to specify regions to be 
                          analyzed; in this case, the op will collapse sample's
                          internal state (see explanation for Data objects)
                          
                          if left None: get the same statistics over the entire
                          sequence; in this case, the op will NOT collapse  
                          sample's internal state
    
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        def token_level_convergence_analysis(data):
            
            recast_data = self._recast(data, where)
            self._prepare_destinations(data, self.dirs.analysis_out)
            
            for sample in recast_data:
                arr = sample.X
                
                if loc is None:
                    freq = self.get_freqs(arr, alphabet)
                    nloc = 'overall'
                    fname = f'{sample.name}_{where}_tokenwise_frequency'
                
                else:                 
                    design = self._infer_design(where)
                    #array internal state has to be collapsed for this calculation
                    if not sample._is_collapsed:
                        msg = f"<frequency_summary> routine will collapse sample {sample.name}'s internal state"
                        self.logger.info(msg)
                        sample._collapse_internal_state()
                    
                    #initialize the frequency array: 3D array to be reduced 
                    #along axis 0 at the end
                    maxlen = self._find_max_len(design, loc)
                    freq = np.zeros((len(design), len(alphabet), maxlen),
                                     dtype=np.float32)
                    
                    for i,template in enumerate(design):                    
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)
    
                        #calculated weighed contributions of each design
                        #to the overall frequency array
                        norm = np.divide(np.sum(row_mask), arr.shape[0])
                        freq[i,:,:len(col_mask)] = norm * np.nan_to_num(
                            
                            self.get_freqs(arr[row_mask][:,col_mask], 
                                           alphabet)
                            
                        )
    
                    #reduce back to a 2D array and plot/save
                    freq = np.sum(freq, axis=0)
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_{where}_reg_{nloc}_tokenwise_frequency'
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.tokenwise_frequency(freq, alphabet, where,
                                                           nloc, basename
                                                          )
                if save_txt:
                    
                    np.savetxt(basename + '.csv',
                               freq,
                               delimiter=',')
                    
            return data
        return token_level_convergence_analysis
    
    
    
    
    
    
    