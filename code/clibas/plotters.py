# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 22:34:18 2021
@author: Alex Vinogradov
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

class SequencingData:
    #Just a container for FastqParser-related plotters
    
    def L_distribution(X, Y, where=None, basename=None):
        
        if not where: where=''
            
        fig = plt.figure(figsize=(18, 6), dpi=300)
        ax = fig.add_subplot(111)
        plt.bar(X, Y, color='#0091b5')
    
        ax.set_ylim(0, 1.02*np.max(Y))
        ax.set_xlim(np.min(X), np.max(X)+1)
        ax.set_xticks(np.linspace(np.min(X), np.max(X)+1, 10))
        ax.set_xticklabels(np.linspace(np.min(X), np.max(X)+1, 10, dtype=int))
        
        ax.set_xlabel('Sequence length', fontsize=30)
        ax.tick_params(axis='both', which='major',  labelsize=25)                                                 
        ax.set_ylabel('Count', fontsize=30)                     
    

        title = f'Distribution of sequence lengths in {where} dataset'
        ax.set_title(title, fontsize=34, y=1.04)
                    
        if basename is not None:                          
            #save png and svg, and close the file
            svg = basename + '.svg'
            png = basename + '.png'
            fig.savefig(svg, bbox_inches = 'tight')
            fig.savefig(png, bbox_inches = 'tight')
            plt.close()    

    def dataset_convergence(C, shannon, where, basename):
    
        y = np.sort(C)
        x = 100 * np.divide(np.arange(y.size), y.size)
    
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)    
        plt.plot(x, y, lw=2.5, c='#3b61b1', antialiased=True)

        ax.set_xlim(0, 101)
        ax.set_xticks(np.arange(0, 125, 25))
        
        ax.set_ylabel(f'{where} sequence count', fontsize=14)
        ax.set_xlabel('Sequence percentile', fontsize=14)
        ax.set_title(f'Sequence-level convergence of {where} dataset', fontsize=16)
        plt.text(x=2, 
                 y=y.max(),
                 s=f'normalized Shannon entropy: {shannon:1.4f}', 
                 size=12,
                 horizontalalignment='left',
                 verticalalignment='center',)
            
        plt.grid(lw=0.5, ls='--', c='slategrey', 
                 dash_capstyle='round', dash_joinstyle='round',
                 antialiased=True, alpha=0.2)  
        
        svg = basename + '.svg'
        png = basename + '.png'
        fig.savefig(svg, bbox_inches = 'tight')
        fig.savefig(png, bbox_inches = 'tight')
        plt.close()      
            
    def conservation(conservation, where, basename):
        
        fig = plt.figure(figsize=(18, 6), dpi=300)
        ax = fig.add_subplot(111)
        plt.plot(conservation, lw=3.5, c='#3b61b1')
    
        y_lim = np.ceil(np.max(conservation))
        ax.set_ylim(0, y_lim)
        ax.set_xlim(0, len(conservation))
        ax.set_xticks(np.linspace(0, len(conservation), 10))
        ax.set_xticklabels(np.linspace(0, len(conservation), 10, dtype=int))
        
        ax.set_xlabel('Sequence index', fontsize=30)
        ax.tick_params(axis='both', which='major',  labelsize=25)                                                 
        ax.set_ylabel('Conservation, bits', fontsize=30)                     
    
        title = f'Token-wise sequence conservation plot for {where} dataset'
        ax.set_title(title, fontsize=34, y=1.04)
                                              
        #save png and svg, and close the file
        svg = basename + '.svg'
        png = basename + '.png'
        fig.savefig(svg, bbox_inches = 'tight')
        fig.savefig(png, bbox_inches = 'tight')
        plt.close()


    def tokenwise_frequency(freq, yticknames, where=None, loc=None, basename=None):
    
        if not where: where=''
        
        if where == 'dna': ylabel = 'Base'
        elif where == 'pep': ylabel = 'Amino acid'
        else: ylabel = 'Token'
            
        figsize = (1 + freq.shape[1] / 2, freq.shape[0] / 2)
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300)
    
        norm = mpl.colors.Normalize(vmin=0, vmax=np.max(freq))
        c = ax.pcolormesh(freq, cmap=plt.cm.Blues, norm=norm, edgecolors='w', linewidths=4)
        cbar = fig.colorbar(c, ax=ax)
    
        cbar.ax.set_ylabel("frequency", rotation=-90, va="bottom", fontsize=22)
        cbar.ax.tick_params(labelsize=20)    
    
        #set ticks
        ax.set_xticks(np.arange(freq.shape[1])+0.5)
        ax.set_yticks(np.arange(freq.shape[0])+0.5)
        ax.set_xticklabels(np.arange(freq.shape[1])+1)
        ax.set_yticklabels(yticknames)
    
        #set labels
        if loc is not None:
            ax.set_xlabel(f'Position inside library region(s) {loc}', fontsize=21)
            
        ax.set_ylabel(ylabel, fontsize=21)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_title(f'Position-wise frequency map for {where} dataset', fontsize=25)
        
        if basename is not None:
            #save png and svg, and close the file
            svg = basename + '.svg'
            png = basename + '.png'
            fig.savefig(svg, bbox_inches = 'tight')
            fig.savefig(png, bbox_inches = 'tight')
            plt.close()
            
        return

    def Q_score_summary(avg, std, loc, basename):
         
        fig = plt.figure(figsize=(18, 6), dpi=300)
        ax = fig.add_subplot(111)
        plt.plot(avg, lw=4, c='#3b61b1')
        plt.plot(avg+std, lw=1, c='#0091b5')
        plt.plot(avg-std, lw=1, c='#0091b5')
        ax.fill_between(np.arange(len(avg)), avg-std, avg+std, color='#0091b5', alpha=0.15)
    
        ax.set_ylim(0, 53)
        ax.set_xlim(0,len(avg))
        ax.set_xticks(np.linspace(0, len(avg), 10))
        ax.set_xticklabels(np.linspace(0, len(avg), 10, dtype=int))
  
        ax.set_xlabel(f'{loc} region(s) index', fontsize=30)
        ax.tick_params(axis='both', which='major',  labelsize=25)                                               
        ax.set_ylabel('Q, average log score', fontsize=30)                     
    
        title = 'Q-score plot'
        ax.set_title(title, fontsize=34, y=1.04)
                                              
        #save png and svg, and close the file
        svg = basename + '.svg'
        png = basename + '.png'
        fig.savefig(svg, bbox_inches = 'tight')
        fig.savefig(png, bbox_inches = 'tight')
        plt.close()   

class Analysis:
    
    def tSNE(embedding, sizes, labels, basename):
        
        import seaborn as sns

        fig = plt.figure(figsize=(8, 8), dpi=300)
        ax = fig.add_subplot(111)
        plt.scatter(embedding[:,0][::-1], embedding[:,1][::-1], 
                    alpha=0.8, 
                    c=labels[::-1],
                    cmap = sns.color_palette("mako", labels.max(), as_cmap=True),
                    marker='o', edgecolors='none', 
                    s=sizes[::-1])
    
        for cluster in labels:
            if not cluster == -1:
                x_coord = np.average(embedding[:,0][labels == cluster])
                y_coord = np.average(embedding[:,1][labels == cluster])
                plt.text(x=x_coord, 
                         y=y_coord,
                         s=f'{cluster+1}', 
                         size=16,
                         weight='bold',
                         alpha=0.8,
                         color='#c38928',
                         horizontalalignment='center',
                         verticalalignment='center')            
    
        ax.set_xlabel('tSNE1, scaled', fontsize=20)
        ax.set_ylabel('tSNE2, scaled', fontsize=20)      
        ax.tick_params(axis='both', which='major',  labelsize=16)                                               
                   
        title = 'tSNE clustering'
        ax.set_title(title, fontsize=24, y=1.04)
                                              
        #save png and svg, and close the file
        svg = basename + '.svg'
        png = basename + '.png'
        fig.savefig(svg, bbox_inches = 'tight')
        fig.savefig(png, bbox_inches = 'tight')
        plt.close()   

    def UMAP_HDBSCAN(Y,
                     labels,
                     C=None, 
                     sample_name=None, 
                     basename=None, 
                     show_annotations=False):
        
        
        cmap = sns.color_palette("husl", labels.max(), as_cmap=True)
        fig = plt.figure(figsize=(8, 8), dpi=300)
        ax = fig.add_subplot(111)
        
        if C is not None:
            sizes = 5000 * np.power(np.divide(C, C.sum()), 0.55)
        else:
            sizes = 5000 * np.power(np.divide(1, Y.shape[0]), 0.55)
        
        if not sample_name: sample_name = 'unnamed sample'
        
        noise = labels == -1
        plt.scatter(Y[:,0][noise][::-1], 
                    Y[:,1][noise][::-1], 
                    alpha=0.6, 
                    c='#070D0D',
                    marker='o', edgecolors='none', 
                    s=sizes[noise][::-1]
                    )     
        
        aint_noise = labels != -1
        plt.scatter(Y[:,0][aint_noise][::-1], 
                    Y[:,1][aint_noise][::-1], 
                    alpha=0.6, 
                    c=labels[aint_noise][::-1],
                    cmap = cmap,
                    marker='o', edgecolors='none', 
                    s=sizes[aint_noise][::-1]
                    )     
      
        if show_annotations:
            for cluster in labels:
                if not cluster == -1:
                    x_coord = np.average(Y[:,0][labels == cluster])
                    y_coord = np.average(Y[:,1][labels == cluster])
                    plt.text(x=x_coord, 
                              y=y_coord,
                              s=f'{cluster}', 
                              size=15,
                              weight='bold',
                              alpha=0.3,
                              color='#323232',
                              horizontalalignment='center',
                              verticalalignment='center')      
         
        plt.axis('off')

        title = f'umap/hdbscan: {sample_name}'
        ax.set_title(title, fontsize=20, y=1.04)
    
        if basename is not None:
            #save png and svg, and close the file
            svg = basename + '.svg'
            png = basename + '.png'
            fig.savefig(svg, bbox_inches = 'tight')
            fig.savefig(png, bbox_inches = 'tight')
            plt.close()
            
        plt.ion
        return

    def ClusteringHyperParams(min_clusters, min_samples, scores, 
                              sample_name=None, 
                              basename=None):
        
        fig = plt.figure(figsize=(6, 7), dpi=300)
        ax = fig.add_subplot(111)
        c = plt.pcolor(scores, cmap=sns.color_palette('mako', as_cmap=True))
        fig.colorbar(c, ax=ax)
        
        ax.set_xlabel('min_sample', fontsize=16)
        ax.set_ylabel('min_cluster', fontsize=16)      
    
        ax.set_xticks(np.arange(len(min_samples)) + 0.5)
        ax.set_yticks(np.arange(len(min_clusters)) + 0.5)
        
        ax.set_xticklabels(min_samples)
        ax.set_yticklabels(min_clusters)    
        
        ax.tick_params(axis='both', which='major',  labelsize=14)                                               
                   
        title = f'hdbscan clustering scores: {sample_name}'
        ax.set_title(title, fontsize=18, y=1.02)    
        
        if basename is not None:
            #save png and svg, and close the file
            svg = basename + '.svg'
            png = basename + '.png'
            fig.savefig(svg, bbox_inches = 'tight')
            fig.savefig(png, bbox_inches = 'tight')
            plt.close()

        return
















