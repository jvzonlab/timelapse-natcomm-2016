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

	worms = glob.glob(path+'*.pickle')
	worms.sort()
	# print(worms)

	divTimes = pd.DataFrame( { } )
	# print(divTimes)

	motherCells = ['a.','b.','c.','1.','2.','3.','4.','5.','6.','t.']

	for worm in worms:
		print( worm )
		df = pickle.load(open(worm,'rb'))

		for side in ['L','R']:
			cdata = df.ix[ (df.rowtype=='cell')&(df.cside==side), ['cname','times','cside'] ]
			finalSeamCells = cdata.ix[cdata.times==np.max(cdata.times.values),'cname'].values
			# print(finalSeamCells)
			cdata.ix[cdata.cname.str.len()==2, 'divCycle'] = 0

			cells = list( set( [ i for i in cdata.cname ] ) )
			cells = [ i for i in cells if len(i) > 2 ]
			cells.sort(key=lambda x: (x[0].isdigit(), x))

			# print(cells)

			# print(cdata)

			for cell in cells:
				# print(cell,( cell in finalSeamCells ))
				if isDaughter(cell, cdata) or ( cell in finalSeamCells ):
					tDiv = findFirstTp( cell, cdata )
					tPrevDiv = findFirstTp( cell[:-1], cdata )

					if len(cell) == 3:
						divCycle = 1
						cdata.ix[ cdata.cname == cell, 'divCycle' ] = divCycle

						divTimes = divTimes.append(pd.Series([worm,side,cell[:2],tDiv,divCycle], index = ['worm','cside','motherCell','tDiv','divCycle']),ignore_index=True)

					elif (tDiv - tPrevDiv) > 4:
						divCycle = cdata.ix[ cdata.cname == cell[:-1], 'divCycle' ].values[0] + 1
						cdata.ix[ cdata.cname == cell, 'divCycle' ] = divCycle
						# print([worm,side,cell,cell[:2],tDiv,divCycle])

						divTimes = divTimes.append(pd.Series([worm,side,cell[:2],tDiv,divCycle], index = ['worm','cside','motherCell','tDiv','divCycle']),ignore_index=True)

					else:
						divCycle = cdata.ix[ cdata.cname == cell[:-1], 'divCycle' ].values[0]
						cdata.ix[ cdata.cname == cell, 'divCycle' ] = divCycle


		
		# find the deviation from the cluster average in the worm
		wormFilter = ( divTimes.worm == worm )
		for cluster in list( set( divTimes.ix[ wormFilter, 'divCycle' ].values ) ):

			clusterFilter = ( divTimes.divCycle == cluster )
			
			mean = np.mean( divTimes.ix[ wormFilter & clusterFilter, 'tDiv' ].values )
			
			divTimes.ix[ wormFilter & clusterFilter, 'deviation' ] = divTimes.ix[ wormFilter & clusterFilter, 'tDiv' ] - mean
			# print(cells)

	# print(divTimes)

	# plot results
	fig, axarr = plt.subplots(int(np.max(divTimes.divCycle)),1,figsize=(5.8,5.8))
	fig.subplots_adjust(left=0.12, right=.95, top=.95, bottom=0.05, hspace = 0.1)

	### plot dotplot
	for cluster in list(set(divTimes.divCycle.values)):

		color = 'black'

		clusterFilter = ( divTimes.divCycle == cluster )

		for idx, mc in enumerate( motherCells ):
			mcFilter = ( divTimes.motherCell.str.startswith( mc ) )
			data = divTimes.ix[ clusterFilter & mcFilter ]

			data.ix[:,'deviation'] = np.floor(data.deviation*10000)/10000
		
			timeLimit = 3
			hist, bins = np.histogram( data.deviation.values, 2*timeLimit*12+1, range = (-timeLimit,timeLimit) )

			bins = bins[:-1]+(bins[1]-bins[0])/2
			step = .06

			print(len(hist),len(bins))

			for val in zip( hist, bins ):
				pos = (np.arange(val[0])-val[0]/2+.5)*step
				for i in np.arange(val[0]):
					axarr[cluster-1].plot( idx-.5+pos[i], val[1], marker = 'o', ms = 3, color = color, mec='none' )

			axarr[cluster-1].plot(idx-.5, np.mean(data.ix[:,'deviation'].values),marker='o',color = 'red',ms=5,mec='none')



	axarr[0].xaxis.tick_top()
	for ax in axarr:
		ax.plot([0,9],[0,0],'--k',lw=.5)
		ax.set_ylim(-2.5,2.5)
		ax.get_yaxis().set_ticks([-2,-1,0,1,2])
		ax.set_xlim(0,9)
		ax.get_xaxis().set_ticks([.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
		ax.invert_yaxis()

	for ax in axarr:
	    for tick in ax.xaxis.get_major_ticks():
	        tick.label.set_fontsize(15)
	    		
	    for tick in ax.yaxis.get_major_ticks():
	        tick.label.set_fontsize(15)

	for ax in axarr:
		ax.xaxis.set_ticklabels([])
	axarr[0].xaxis.set_ticklabels([.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
	realSeamCellNames = ['H1','H2','V1','V2','V3','V4','V5','V6','T']
	axarr[0].set_xticklabels( realSeamCellNames, fontsize = 15 )

	axarr[0].set_ylabel('T_1-<T_1>',fontsize = 15)
	axarr[1].set_ylabel('T_2-<T_2>',fontsize = 15)
	axarr[2].set_ylabel('T_3-<T_3>',fontsize = 15)
	axarr[3].set_ylabel('T_4-<T_4>',fontsize = 15)	






if __name__ == '__main__':

    # path='Y:\\Images\\150617_JVZ20_wrt2GFP\\'
    path = 'Y:\\Images\\150907_JR667_250x250x20\\'
    # path = 'Y:\\Images\\EMS3_alldata\\'
    
    plotDivisions(path)

    plt.show()