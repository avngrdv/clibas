# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 19:20:25 2022
@author: Alex Vinogradov

Module contains holoviews plotting instruments (as opposed to the classic 
matplotlib ones in plotters.py).

At some stage the two should be merged, and everything converted to hv, 
but for now hv is kept separate for cleanliness. 
"""

import numpy as np
import holoviews as hv
from bokeh.models import HoverTool
import seaborn as sns
hv.extension('bokeh')

def _get_heatmap(X, labels, cluster=None, alphabet=None):
    
    from clibas.misc import get_freqs
    
    if cluster is not None:
        freq = get_freqs(X[labels == cluster], alphabet=alphabet)
    else:
        freq = get_freqs(X, alphabet=alphabet)
        
    hmap = [(i, j, freq[j, i]) 
            for i in np.arange(freq.shape[1]) 
            for j in np.arange(freq.shape[0])
           ]
    
    return hmap

def _get_palette(clusters, bw=False):
    
    n_clusters = np.unique(clusters).size
    palette = sns.color_palette("husl", n_clusters).as_hex()
    
    if 0 in clusters:
        palette[0] = '#323232'
        
    if bw:
        palette = ['#323232' for c in palette]
        
    return palette

def _XYC_unpack(d):
    umap1 = d['Y'][:,0]
    umap2 = d['Y'][:,1]
    seqs = np.array([''.join(x) for x in d['X']])
    sizes = 80 * np.power(np.divide(d['C'], d['C'].sum()), 0.25)
    return umap1, umap2, seqs, sizes

def hv_umap_embedding(seqs=None, 
                      umap1=None, 
                      umap2=None, 
                      C=None, 
                      sizes=None, 
                      labels=None,
                      clusters=None,
                      lims=None,
                      sname='',
                      bw=False): #bw: black and white
    
    d = {'x': umap1,
         'y': umap2,
         'size': sizes,
         'count': C,
         'cluster': labels,
         'seqs': seqs,
         'top': np.arange(seqs.size) + 1
        }
    
    #make a custom palette:
    #no idea how to do it with bokeh using the husl colorwheel
    palette = _get_palette(clusters, bw=bw)
    
    #formatting for hover tooltips
    TOOLTIPS = """
        <div>
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold; ">@top: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@seqs</span>
            </div>
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold;">Count: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@count</span>
            </div>  
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold;">Cluster: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@cluster</span>
            </div>            
        </div>
    """
    
    if lims is not None:
        xlim = lims[0]
        ylim = lims[1]
    else:
        xlim = (umap1.min() - 0.25, umap1.max() + 0.25)
        ylim = (umap2.min() - 0.25, umap2.max() + 0.25)
    
    title = f'{sname} UMAP embeddings'
    
    scatter = hv.Scatter(d, 
                         vdims=['y', 'size', 'cluster', 'seqs', 'count', 'top'],
                         label=title
                        )
    
    scatter = scatter.opts(height=600,
                           width=600,
                           yaxis=None, 
                           xaxis=None,
                           xlim=xlim,
                           ylim=ylim
                          )
    
    scatter = scatter.opts(size='size', 
                           alpha=0.7, 
                           line_width=0, 
                           color='cluster',
                           cmap=palette,
                           default_tools = []
                          )
    
    scatter = scatter.opts(tools=[
                                  HoverTool(tooltips=TOOLTIPS),
                                  'reset',
                                  'save',
                                  'box_zoom',
                                  'wheel_zoom'
                                 ], 
                           nonselection_fill_alpha=0.1
                          )
    return scatter

def hv_cluster_bars(cluster_summary, clusters):
    
    d = {
         'size': cluster_summary['Cluster size'],
         'cluster': cluster_summary['Cluster number'],
         'purity': cluster_summary['Cluster purity'],
         'score': cluster_summary['Cluster score']
        }

    def hook(plot, element):
        plot.handles['xaxis'].major_label_text_font_size = '0pt'
        plot.handles['xaxis'].axis_label_text_font_size = '14pt'
        plot.handles['yaxis'].axis_label_text_font_size = '14pt'
        
    #make a custom palette:
    #no idea how to do it with bokeh using the husl colorwheel
    palettte = _get_palette(clusters, bw=False)

    TOOLTIPS = """
        <div>
          <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold; ">Cluster: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@cluster</span>
            </div>       
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold; ">Score: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@score</span>
            </div>
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold;">Purity: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@purity</span>
            </div>  
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold;">Size: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@size</span>
            </div>            
        </div>
    """

    bars = hv.Bars(d, 
                   hv.Dimension('cluster'),
                   vdims=['score', 'purity', 'size'],
                   label='Clustering quality'
                  )
    
    bars = bars.opts(height=200,
                     width=600,
                     ylim=(0, cluster_summary['Cluster score'].max() + 0.5), 
                     cmap=palettte, 
                     color='cluster',
                     line_width=0
                    )
    
    bars = bars.opts(tools=[HoverTool(tooltips=TOOLTIPS), 
                            'reset',
                            'save',
                            'box_zoom',
                            'wheel_zoom'
                           ],
                     nonselection_fill_alpha=0.1
                    )
    
    bars = bars.redim.values(**{'cluster': clusters})

    bars = bars.opts(ylabel='Cluster score',
                     xlabel='Cluster (hover for details)',
                     hooks=[hook],
                    )
    return bars

def hv_cluster_freqs(freq_heat_map, alphabet=None, x_dim_size=None):
    
    hm = hv.HeatMap(freq_heat_map, vdims=['z'])

    xticks = [(int(i), i+1) for i in np.arange(x_dim_size)]
    yticks = [(i, token) for i,token in enumerate(alphabet)]

    cbar_opts = {'height': 200,
                 'width': 15, 
                 'location': 'top_right'
                }

    TOOLTIPS = """
        <div>
            <div>
                <span style="font-size: 14px; color: #905c54; font-family: Arial; font-weight: bold; ">freq: </span>
                <span style="font-size: 15px; color: #323232; font-family: Consolas; font-weight: bold;">@z{(0.00)}</span>
            </div>          
        </div>
    """
    
    hm = hm.opts(xticks=xticks, 
                 yticks=yticks,
                 xlabel='Position',
                 ylabel='Token',
                 title='Sequence conservation'
                )
    
    hm = hm.opts(height=600,
                 width=35 * x_dim_size, 
                 colorbar=True, 
                 colorbar_opts=cbar_opts,
                 fontsize={'title': 9}
                )
    
    hm = hm.opts(tools=[HoverTool(tooltips=TOOLTIPS),
                        'save',
                        'reset',           
                        'box_zoom',
                        'wheel_zoom'
                       ],
                 fontscale=1.4,
                )
    return hm
	

def hdbumap_holomap_triplet(d,
                            alphabet=None):

    umap1, umap2, seqs, sizes = _XYC_unpack(d)
    
    cluster_summary = d['cluster_summary']
    labels = d['labels']
    
    clusters = np.array(cluster_summary['Cluster number'])
    x_dim_size = d['X'].shape[-1]

    scs_lims = [
                (umap1.min() - 0.25, umap1.max() + 0.25),
                (umap2.min() - 0.25, umap2.max() + 0.25)
               ]
    
    #initialize the objects
    hms = dict()
    brs = dict()
    scs = dict()

    #first entry is all clusters together, labelled as -7
    hmap = _get_heatmap(d['X'], labels, alphabet=alphabet)
    
    hms.update({'All together': hv_cluster_freqs(hmap,
                                                 alphabet=alphabet,
                                                 x_dim_size=x_dim_size
                                                 )
               })
       
    brs.update({'All together': hv_cluster_bars(cluster_summary, clusters)})
    scs.update({'All together': hv_umap_embedding(seqs=seqs, 
                                                  umap1=umap1,
                                                  umap2=umap2, 
                                                  C=d['C'], 
                                                  sizes=sizes,
                                                  labels=labels,
                                                  clusters=clusters,
                                                  sname=d['name'],
                                                  lims=scs_lims
                                                 )
               })

    for i in clusters:
        max_dig = len(str(clusters.max()))
        entry = f'cluster {i:0{max_dig}}'

        hms.update({entry: hv_cluster_freqs(_get_heatmap(d['X'], 
                                                         labels,
                                                         i, 
                                                         alphabet=alphabet
                                                        ),
                                        
                                            alphabet=alphabet, 
                                            x_dim_size=x_dim_size
                                           ) 
                   })
    
        brs.update({entry: hv_cluster_bars(cluster_summary[
                                           cluster_summary['Cluster number'] == i
                                                          ],
                                           clusters
                                          ) 
                  })
        
        scs.update({entry: hv_umap_embedding(
                                             seqs=seqs[labels == i], 
                                             umap1=umap1[labels == i],
                                             umap2=umap2[labels == i], 
                                             C=d['C'][labels == i], 
                                             sizes=sizes[labels == i],
                                             labels=labels[labels == i],
                                             clusters=clusters,
                                             lims=scs_lims
                                            )
                  })

    return hms, brs, scs

def hdbumap_analysis_dashboard(d,
                               alphabet=None,
                               fname=None):
                            
    hms, brs, scs = hdbumap_holomap_triplet(d, alphabet=alphabet)
    
    a = hv.HoloMap(hms, kdims='Cluster number').opts(default_tools = [])
    b = hv.HoloMap(brs, kdims='Cluster number').opts(default_tools = [])
    c = hv.HoloMap(scs, kdims='Cluster number').opts(default_tools = [])
    
    L = c + a + b
    L.cols(2)
    
    hv.output(widget_location='bottom')
    L.opts(shared_axes=False)
    
    fname += '.html'
    hv.save(L, fname)
    return

def single_manifold_embedding_dashboard(tup, fname=None):
    
    scs = list()
    
    #first pass to determine the common coordinate grid        
    x_min = [d['Y'][:,0].min() for d in tup]
    x_max = [d['Y'][:,0].max() for d in tup]
    y_min = [d['Y'][:,1].min() for d in tup]
    y_max = [d['Y'][:,1].max() for d in tup]
    
    scs_lims = (
                (min(x_min) - 0.25, max(x_max) + 0.25),
                (min(y_min) - 0.25, max(y_max) + 0.25),
               )     

    #now unpack and prepare plots
    for d in tup:
    
        umap1, umap2, seqs, sizes = _XYC_unpack(d)
        clusters = np.array(d['cluster_summary']['Cluster number'])   
        
        scatter  =   hv_umap_embedding(seqs=seqs, 
                                       umap1=umap1, 
                                       umap2=umap2, 
                                       C=d['C'], 
                                       sizes=sizes, 
                                       labels=d['labels'],
                                       clusters=clusters,
                                       lims=scs_lims,
                                       sname=d['name'],
                                       bw=True
                                      )
        scs.append(scatter)

    L = hv.Layout(scs).cols(2)
    fname += '.html'
    hv.save(L, fname)
    return

