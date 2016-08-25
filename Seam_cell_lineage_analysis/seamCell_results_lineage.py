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

def drawLineage( path, worm, maxTp = 100 ):

    def findCellData( cname, cdata, ax ):
        cell = { 'cname': cname, 'firsttp': np.nan, 'lasttp': np.nan, 'isDaughter': False }

        ft = cdata.ix[ cdata.cname == cname, 'times' ].values

        # no daughter, but granddaughter present
        if (not isDaughter(cname, cdata)) & isLineageContinuing(cname, cdata):
        	# find the next cells in the lineage
            i=1
            foundnewmcs = True
            while foundnewmcs:
                newmcs = list( set( cdata.ix[ cdata.cname.str.startswith(cname) & ( cdata.cname.str.len() == (len(cname)+i) ), 'cname'] ) )
                newmcs.sort()
                i+=1
                foundnewmcs = len(newmcs)==0
            # start a new mothercell lineage
            for newmc in newmcs:
            	drawSingleLineage(newmc, ax, cdata, finalSeamCells)

        # if no time points found, look at the sister
        elif len( ft ) == 0:
            if cname[-1] == 'p':
                sisterCellName = cname[:-1] + 'a'
            else:
                sisterCellName = cname[:-1] + 'p'
            sisterCell = findCellData( sisterCellName, cdata, ax )
            cell[ 'firsttp' ] = sisterCell[ 'firsttp' ]
            cell[ 'lasttp' ] = cell[ 'firsttp' ] + 3
            return cell
        
        # the first tp should be 0 for original mothercells
        if len(cname) == 2:
            cell[ 'firsttp' ] = 0
        else:
            cell[ 'firsttp' ] = np.min( ft )

        # set if a daughter exists
        cell['isDaughter'] = isDaughter( cname, cdata )

        # if a duaghter exists, set her first tp as last tp of mother
        if cell['isDaughter']:
            lt = cdata.ix[ ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() > len(cname) ), 'times' ].values
            lt = np.min( lt )
        # if this is a final seam cell, extend the line
        elif cell['cname'] in finalSeamCells:
            lt = maxTp
        # if it's not, make the line shorter
        else:
            lt = np.min( ft )+3
        cell[ 'lasttp' ] = lt
        return cell

    def isDaughter(cname, cdata):
        return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() == (len(cname)+1) ) ).any()

    def isLineageContinuing(cname, cdata):
        return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() > len(cname) ) ).any()

    def drawHorLine(cell,pos,ax,cdata, finalSeamCells):
        mean = np.mean(pos)
        color = ['gray','gray']
        if isDaughter( cell['cname']+'a', cdata ) or (cell['cname']+'a') in finalSeamCells:
            color[0] = 'black'
        if isDaughter( cell['cname']+'p', cdata ) or (cell['cname']+'p') in finalSeamCells:
            color[1] = 'black'
       	ax.plot( [ pos[0], mean ], [ cell['lasttp'], cell['lasttp'] ], '-', color = color[0], lw=1.5 )
        ax.plot( [ mean, pos[1] ], [ cell['lasttp'], cell['lasttp'] ], '-', color = color[1], lw=1.5 )

    def drawSingleLineage(cname,ax,cdata,finalSeamCells):
        cell = findCellData( cname, cdata, ax )
        # print('found data for cell:', cname, cdata.cside.values[0])
        pos = positions[ cell['cname'][:2] ]
        n = 1.
        for ap in cell['cname'][2:]:
        	pos += 0.4 * ( ( ap=='p' ) - ( ap=='a' ) ) / (n)
        	n+=1.
        if cell['isDaughter'] or ( cell['cname'] in finalSeamCells ):
            color = 'black'
        else:
            color = 'gray'
        ax.plot( [ pos, pos ], [ cell['firsttp'], cell['lasttp'] ], '-',color = color, lw=1.5 )
        if cell['isDaughter']:
            drawHorLine(cell,[pos-0.4/n,pos+0.4/n],ax, cdata, finalSeamCells)
            for c in [ cname+'a', cname+'p' ]:
                drawSingleLineage( c, ax, cdata, finalSeamCells)
        else:
            return


    # load the data    
    if os.path.exists( os.getcwd()+'\\worm'+worm+'.pickle' ):
        loadPath = os.getcwd()+'\\worm'+worm+'.pickle'
    else:
        loadPath = path + 'worm' + worm + '.pickle'
    dfOriginal = pickle.load(open(loadPath,'rb'))
    df = dfOriginal.ix[ dfOriginal.times < maxTp ]



    for side in ['L','R']:
        # # setup the figure
        fig = plt.figure(figsize = (5.8,3.5))
        ax = fig.add_subplot(111)
        ax.set_xlim(0.5,11.5)
        ax.set_ylim(0,maxTp)
        ax.invert_yaxis()
        fig.subplots_adjust(left=0.12, right=.99, top=.90, bottom=0.05)

        # format plot
        for tl in ax.get_yticklabels():
            tl.set_fontsize(15) 
        ax.get_xaxis().set_visible(False)

        ax.text(-.5,-1,side,fontsize=15)

        # draw lines for ecdysis

        ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )

        index = np.where( ecd[:,0] == float(worm[1:]) )
        mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
        print(mintp)
        lethtp = ecd[index, 2:][0][0] - mintp - 1
        print(lethtp)
        # print(ecd[index,1:][0][0])

        for tidx in lethtp:
            if ( tidx > 0 ):
                if ( dfOriginal.ix[(dfOriginal.rowtype=='body')&(dfOriginal.tidx==tidx),'times'].values[0] < maxTp ):
                    tp = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]
                    ax.plot([.5,11.5],[tp,tp],'--',color=[.7,.7,.7])

        # filter only the data of cells
        cdata = df.ix[ (df.rowtype=='cell')&(df.cside==side), ['cname','times','cside'] ]
        finalSeamCells = cdata.ix[cdata.times==np.max(cdata.times.values),'cname'].values

        # motherCells = list( set( cdata.ix[ cdata.cname.str.len()==2, 'cname'] ) )
        # motherCells.sort(key=lambda x: (x[0].isdigit(), x))
        motherCells=['a.','b.','c.','1.','2.','3.','4.','5.','6.','t.']
        print(motherCells)

        positions = { i[0]: i[1] for i in zip(motherCells,np.append( np.arange(1,len(motherCells)), len(motherCells)+.5 ) ) }
        realNames = { i[0]: i[1] for i in zip(motherCells,['H0','H1','H2','V1','V2','V3','V4','V5','V6','T'])}

        for mc in motherCells:
            ax.text(positions[mc]-.25,-1,realNames[mc],fontsize=15)
            drawSingleLineage(mc, ax, cdata, finalSeamCells)



if __name__ == '__main__':

    path='Y:\\Images\\150907_JR667_250x250x20\\'
    worms = ['C19']

    path = 'X:\\Michael\\160317_EMS7\\'
    # worms = ['C15','C16','C19','C20']
    worms = ['C16']

    # path = 'Y:\\Images\\151111_EMS3_250x250x20\\'
    # worms = ['C03','C07','C12']

    # path = 'Y:\\Images\\150824_JVZ20_wrt2GFP_250x250x20\\'
    # worms = ['C01']

    for worm in worms:    
        drawLineage(path,worm,maxTp=50)

    plt.show()