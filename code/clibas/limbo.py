# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:58:01 2022

@author: a_vin
"""

    def tSNE_analysis(self, where=None, top_n=1000, cluster_fasta=False):
        '''
        For each sample in Data, compute tSNE embeddings for the dataset 
        specified by 'where'. Cluster the results (HDBSCAN) and summarize. 
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on
    
                   top_n: the number of top entries (by count) to analyze   
    
           cluster_fasta: True/False. if True, create a .fasta file for
                          sequences comprising each cluster. 
          						  						   
        Returns:
                Data object (no transformation)
        '''        
    
        self._where_check(where)
        if not isinstance(top_n, int):
            msg = f'<tSNE_analysis> routine expected to receive parameter top_n as as int; received: {type(top_n)}'
            self.logger.error(msg)  
            raise ValueError(msg)        
    
        if not isinstance(cluster_fasta, bool):
            msg = f'<tSNE_analysis> routine expected to receive parameter cluster_fasta as as bool; received: {type(cluster_fasta)}'
            self.logger.error(msg)  
            raise ValueError(msg)        
    
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            import hdbscan
        except:
            msg = 'Failed to import the libraries necessary for <tSNE_analysis> op. . .'
            self.logger.error(msg)
            raise ImportError(msg)
    
        def _cluster_writer(df, destination, sname):
            clusters = np.unique(df['cluster'])
            for cluster in clusters:
                fname = f'{sname}_{where}_HDBSCAN_cluster_{cluster}.fasta'
                path = os.path.join(destination, fname)
                
                with open(path, 'w') as f:
                    for entry in df[df['cluster'] == cluster].iterrows():
                        f.write(f'>cluster_{cluster}_count_{entry[1]["counts"]}\n')  
                        f.write(f'{entry[1]["sequence"]}\n')
            return        
    
        def tSNE_embedding(data):
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                #get an array of top_n entries, X; C - counts array
                X, C = np.unique(arr, return_counts=True, axis=0)
                ind = np.argsort(C)[::-1][:top_n]
                X = X[ind]
                C = C[ind]
                
                #get one-hot encodings of the datasets
                tokens = np.unique(arr)
                nX = np.vectorize(lambda x: np.where(tokens == x)[0][0])(X)
                nX = nX.astype(np.int8).ravel()
                
                fX = np.zeros((X.size, tokens.size))
                fX[np.arange(nX.size), nX] = 1
                fX = np.reshape(fX, (X.shape[0], -1))
                
                #embed the representations and scale
                embedder = TSNE(n_components=2,
                                perplexity=30, 
                                learning_rate=200,
                                init='random')
                
                Y = embedder.fit_transform(fX.astype(np.float32))
                Y = StandardScaler().fit_transform(Y)
                
                #run HDBSCAN clustering
                hdb = hdbscan.HDBSCAN(cluster_selection_epsilon=0.02,
                                      cluster_selection_method = 'eom')
                hdb.fit(Y)
                
                #plot the results
                destination = os.path.join(self.dirs.parser_out, sample.name, 'tSNE_analysis')
                if not os.path.isdir(destination):
                    os.makedirs(destination)       
                    
                fname = f'{sample.name}_{where}_tSNE'
                basename = os.path.join(destination, fname)         
                sizes = 8000 * np.power(np.divide(C, arr.shape[0]), 0.6)
                Plotter.Analysis.tSNE(Y, sizes, hdb.labels_, basename)
                
                #dump analysis into a .csv
                entries = np.array([''.join(x) for x in X])
                clusters = hdb.labels_ + 1
                
                d = {'sequence': entries,
                     'counts': C,
                     'cluster': clusters}
                
                df = pd.DataFrame.from_dict(d)
                df = df.sort_values(['cluster', 'counts'], ascending=False)
                fname = f'{sample.name}_{where}_HDBSCAN_clustering.csv'
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=',', index=False)
                
                if cluster_fasta:
                    _cluster_writer(df, destination, sample.name)
                
            return data
        return tSNE_embedding
    