# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 15:43:57 2015

@author: gritti
"""

import glob
from timelapseFun import *
from tifffile import *
import numpy as np
from PIL import Image
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['pdf.fonttype'] = 42

# import seaborn as sns
# sns.set_style('ticks')


def plotDivisions(path):

    def isDaughter( cname, cdata ):
        return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() == (len(cname)+1) ) ).any()

    def findFirstTp( cname, cdata ):
        # if a duaghter exists, set her first tp as last tp of mother
        lt = cdata.ix[ ( cdata.cname == cname ), 'times' ].values
        lt = np.min( lt )
        return lt

    def findLastTp( cname, cdata ):
        # if a duaghter exists, set her first tp as last tp of mother
        lt = cdata.ix[ ( cdata.cname == cname ), 'times' ].values
        lt = np.max( lt )
        return lt

    worms = glob.glob(path+'*.pickle')
    worms.sort()
    # print(worms)
    motherCells = ['b.','c.','1.','2.','3.','4.','5.','6.','t.']

    rows = []
    for i in worms:
        rows.append(i[-10:-7]+'_L')
        rows.append(i[-10:-7]+'_R')
    print(rows)
    divTimes = pd.DataFrame( { 'worm': [ i for i in rows ] } )
    for i in motherCells:
        divTimes.ix[:,i] = np.nan
    print(divTimes)

    for worm in worms:
        print( worm )
        df = pickle.load(open(worm,'rb'))

        for side in ['L','R']:
            cdata = df.ix[ (df.rowtype=='cell')&(df.cside==side), ['cname','times','cside'] ]

            for cell in motherCells:
                # print(cell,( cell in finalSeamCells ))
                tDiv = findFirstTp( cell+'p', cdata )

                divTimes.ix[ ( divTimes.worm == worm[-10:-7]+'_'+side ), cell ] = tDiv

    divTimes.ix[ 16, 'b.': ] = np.mean( divTimes.ix[ :, 'b.': ], 0 )
    divTimes.ix[ 16, 'worm' ] = 'average'

    print(divTimes)

    # assign order sequence
    divSequence = pd.DataFrame( { 'worm': [ i for i in rows ] } )
    for i in motherCells:
        divSequence.ix[:,i] = np.nan
    divSequence.ix[ 16, 'worm' ] = 'average'

    for idx, data in divTimes.iterrows():
        wormname = data.worm

        tdivdata = data['b.':]
        tdivdata.sort()
        print('\n',wormname)

        sequence = np.ones((1,9))
        for t in tdivdata:
            sequence += data['b.':]>t
        divSequence.ix[ divSequence.worm == wormname, 'b.': ] = sequence

    print(divSequence)

    fig = plt.figure(figsize = (5.8,3.5))
    ax = fig.add_subplot(111)
    ax.set_xticks(np.arange(0,17))
    ax.set_yticks(np.arange(0,10))
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    fig.subplots_adjust(left=0.12, right=.85, top=.95, bottom=0.05)
    img = ax.imshow(divSequence.ix[ :, 'b.': ].T, interpolation = 'nearest', cmap = cm.get_cmap('gist_rainbow',9), vmin = 1, vmax = 9)
    # plt.clim(0,18)
    ax1 = fig.add_axes([.9,.15,.05,.7])
    fig.colorbar(img,cax=ax1)

if __name__ == '__main__':

    # path='Y:\\Images\\150617_JVZ20_wrt2GFP\\'
    path = 'Y:\\Images\\150907_JR667_250x250x20\\'
    
    plotDivisions(path)

    plt.show()
