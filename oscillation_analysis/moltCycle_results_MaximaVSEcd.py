# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 17:12:50 2015

@author: gritti
"""

import glob
import numpy as np
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
import matplotlib as mpl
from matplotlib import cm
from scipy.ndimage import filters
from scipy.signal import argrelextrema
import pandas as pd
mpl.rcParams['pdf.fonttype'] = 42

path = 'Y:\\Images\\150914_mlt10_250x250x20\\'
worms = ['C01','C03','C04','C05','C07']

color = np.array( [cm.jet(int(i))[:3] for i in np.arange(len(worms))*254./(len(worms)-1)] )

# setup the figure
fig = plt.figure(figsize = (5.8,5.8))
ax = fig.add_subplot(111)
tMax = 18
ax.set_xlim(5,tMax)
ax.set_ylim(5,tMax)
fig.subplots_adjust(left=0.12, right=.97, top=.97, bottom=0.12)

# format plot
for tl in ax.get_yticklabels():
	tl.set_fontsize(18) 
for tl in ax.get_xticklabels():
	tl.set_fontsize(18) 
ax.plot( [0, tMax], [0, tMax], '--k' )


for idx, worm in enumerate( worms ):

	# load worm data
	loadPath = path + 'worm' + worm + '.pickle'
	df = pickle.load(open(loadPath,'rb'))
	df.exp = df.exp / 1000

	# load ecdysis data
	ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )
	index = np.where( ecd[:,0] == float(worm[1:]) )
	lethtidx = ecd[index, 2:][0][0] - np.min(ecd[index, 1:][0][0]) - 1
	lethtp = [ df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0] for tidx in lethtidx ]

	# group data by time and remove missing datapoint
	groups = df.groupby(df.times)
	data = groups.mean().ix[ (groups.mean().tidx >= 0) ]
	data = data.ix[pd.notnull(data.exp)]
	data.ix[:,'times'] = data.index
	data = data.reset_index(drop=True)

	# compute finer 1d interpolation
	times = np.linspace( np.min(data.times), np.max(data.times), 60 * ( np.max(data.times) - np.min(data.times) ) )
	interpolation = ip.interp1d(data.times,data.exp,kind='linear')
	interpexp = interpolation(times)

	# compute gaussian filter of data
	filtdata = pd.DataFrame( {
		'times': times,
		'filtexp': filters.gaussian_filter1d(interpexp,60)
		} )
	# plt.plot(filtdata.times,filtdata.filtexp)

	# find maxima
	maxexptp = filtdata.times[ argrelextrema( filtdata.filtexp.values, np.greater )[0] ].values

	# assign peaks to lethargus
	tpdata = pd.DataFrame( {'tEcd': lethtp, 'tMaxExp': np.nan } )

	for tp in lethtp:
		dist = np.abs( maxexptp - tp )
		closestMaxTp = maxexptp[ np.where( dist == np.min(dist) ) ]
		# print(closestMaxTp)
		tpdata.ix[ tpdata.tEcd == tp, 'tMaxExp' ] = closestMaxTp


	ax.plot( tpdata.tEcd, tpdata.tMaxExp, 'o', 
		fillstyle = 'full', 
		mew=0.1, 
		ms = 8, 
		# mec = 'black', 
		mfc = color[idx], 
		alpha = .8 )

	ax.set_xlabel('Ecdysis', fontsize=18)
	ax.set_ylabel('Peak', fontsize=18)

plt.show()
