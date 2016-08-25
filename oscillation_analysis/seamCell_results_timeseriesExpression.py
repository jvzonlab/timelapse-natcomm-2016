# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:22:21 2015

@author: gritti
"""

import glob
import numpy as np
import os.path
import matplotlib.pyplot as plt
import pickle
import matplotlib as mpl
from matplotlib import cm
import pandas as pd
from scipy.ndimage import filters
import scipy.interpolate as ip
mpl.rcParams['pdf.fonttype'] = 42

def plotExpressionOneWorm( axes, path='Y:\\Images\\150617_JVZ20_wrt2GFP\\', worm='C01', 
    plotCells = True, plotMean = True, plotSmooth = True, plotLeth = True, plotNSeam = True, color = [0,0,0,1] ):
    
    # load the data    
    if os.path.exists( os.getcwd()+'\\worm'+worm+'.pickle' ):
        loadPath = os.getcwd()+'\\worm'+worm+'.pickle'
    else:
        loadPath = path + 'worm' + worm + '.pickle'
    df = pickle.load(open(loadPath,'rb'))
    df.wrt2exp = df.wrt2exp / 1000
        
    # plot the gray dots for single cells
    if plotCells:
        axes.plot(df.times,df.wrt2exp,'o',color = [.8,.8,.8], alpha = .7, mec = None, mew = 0.01, markersize = 4.)
    
    ## PLOT THE MEAN
    # plot the mean expression level binned every 1/N hour
    bins = np.arange(0,50,.5)
    groups = df.groupby(pd.cut(df.times,bins))
    
    if plotMean:
        pl1, = ax.plot(groups.mean().times[pd.notnull(groups.mean().wrt2exp)],
            groups.mean().wrt2exp[pd.notnull(groups.mean().wrt2exp)],'-',lw=2,color = color)

    ## PLOT THE FILTERED DATA

    if plotSmooth:
        # group data by time and remove missing datapoint
        groups = df.groupby(df.times)
        data = groups.mean().ix[ (groups.mean().tidx >= 0) ]
        data = data.ix[pd.notnull(data.wrt2exp)]
        data.ix[:,'times'] = data.index
        data = data.reset_index(drop=True)

        # compute finer 1d interpolation
        times = np.linspace( np.min(data.times), np.max(data.times), 60 * ( np.max(data.times) - np.min(data.times) ) )
        interpolation = ip.interp1d(data.times,data.wrt2exp, kind='linear')
        interpexp = interpolation(times)

        # compute gaussian filter of data
        filtdata = pd.DataFrame( {
            'times': times,
            'filtexp': filters.gaussian_filter1d(interpexp,60)
            } )
        pl1, = ax.plot(filtdata.times,
                filtdata.filtexp,'-',lw=3,color = color)

    ## PLOT NUMBER OF SEAM CELLS
    if plotNSeam:
        # create right axes
        ax2 = axes.twinx()

        times = []
        scn = []
        for t in df.times:
            if len(df.ix[df.times==t]) > 1:
                times.append(t)
                scn.append(len(df.ix[df.times==t]) - 1)
        pl2, = ax2.plot( times, scn, 'b-', lw = 2 )

        for tl in ax2.get_yticklabels():
            tl.set_fontsize(15)    

        ax2.set_ylim(0,30*5)
        ax2.plot([0,45],[30,30],'--',color = pl2.get_color(), lw = 2)
        ax2.set_ylabel('# seam cells', fontsize = 15)
        ax2.set_yticks(np.arange(0,30+1,10))    
        ax2.yaxis.label.set_color(pl2.get_color())
        ax2.tick_params(axis='y', colors=pl2.get_color())

    ## PLOT LETH DATA
    if plotLeth:
        ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )

        index = np.where( ecd[:,0] == float(worm[1:]) )
        lethtp = ecd[index, 2:][0][0] - np.min(ecd[index, 1:][0][0]) - 1

        for tidx in lethtp:
            tp = df.ix[ (df.rowtype=='body') & (df.tidx==tidx), 'times' ].values[0]
            axes.plot([tp,tp],[0,20],'--', color = color, lw=2)



if __name__ == '__main__':

    ### SETUP THE FIGURE
    fig = plt.figure(figsize = (5.8,3.8))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.12, right=.97, top=.97, bottom=0.16)
    for tl in ax.get_xticklabels():
        tl.set_fontsize(18)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(18)     
    ax.set_ylim([0,20])
    ax.set_xlim([0,45])
    ax.set_xlabel('Hours after hatching', fontsize = 18)
    ax.set_ylabel('Nuclear fluorescence [a.u.]', fontsize = 18)

    ### DEFINE ALL THE WORMS
    path = 'Z:\\150824_JVZ20_wrt2GFP_250x250x20\\'
    worms = ['C02','C03','C04','C06','C10','C12','C14']

    color = np.array( [cm.jet(int(i))[:3] for i in np.arange(len(worms))*254./(len(worms)-1)] )

    # ## PLOT WORMS DATA
    # for idx, worm in enumerate( worms ):
    #     print( path, worm )
    #     plotExpressionOneWorm( ax, path, worm, plotCells=False, plotMean=False, plotSmooth=True, plotLeth = True, plotNSeam = False, color = color[idx])

    # PLOT SINGLE WORM DATA
    worm = 'C14'
    print( path, worm )
    plotExpressionOneWorm( ax, path, worm, plotCells=True, plotMean=False, plotSmooth=True, plotLeth = True, plotNSeam = False, color = [0,0,0,1])

    plt.show()
