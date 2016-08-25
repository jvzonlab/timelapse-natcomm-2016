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
from matplotlib.ticker import NullFormatter
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib.ticker import AutoMinorLocator


path = 'Y:\\Images\\150907_JR667_250x250x20\\'
worms = glob.glob(path+'*.pickle')
worms.sort()

### setup figure for the timeseries
fig1 = plt.figure(figsize=(2,2*35.677/42.84))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)


### plot Z diameter

mothercell = '3.'
# average = [[3,9,15,21,27,33],[0,0,0,0,0,0]]
oldaverage = 0
for worm in worms:
	print( worm )
	df = pickle.load(open(worm,'rb'))

	cdata = df.ix[ (df.rowtype=='cell') ]
	tdata = df.ix[ (df.rowtype=='body'), ['times', 'tidx'] ]
	# print(tidx)
	data = [ [], [] ]

	for idx, trow in tdata.iterrows() :
		cellsL = cdata.ix[ cdata.cname.str.startswith(mothercell) & ( cdata.tidx == trow.tidx ) * ( cdata.cside == 'L' ) ]
		cellsR = cdata.ix[ cdata.cname.str.startswith(mothercell) & ( cdata.tidx == trow.tidx ) * ( cdata.cside == 'R' ) ]

		if (len(cellsL)>0) and (len(cellsR)>0) :
			# print(cellsL, cellsR)

			Z1  = cellsL.ix[ cellsL.cXpos == np.max( cellsL.cXpos ), 'cZpos' ].values[0]
			Z2  = cellsR.ix[ cellsR.cXpos == np.max( cellsR.cXpos ), 'cZpos' ].values[0]

			X1  = cellsL.ix[ cellsL.cXpos == np.max( cellsL.cXpos ), 'cXpos' ].values[0]
			X2  = cellsR.ix[ cellsR.cXpos == np.max( cellsR.cXpos ), 'cXpos' ].values[0]

			Y1  = cellsL.ix[ cellsL.cXpos == np.max( cellsL.cXpos ), 'cYpos' ].values[0]
			Y2  = cellsR.ix[ cellsR.cXpos == np.max( cellsR.cXpos ), 'cYpos' ].values[0]

			# data[0].append( np.sqrt( ( (Z1-Z2)*2.5 )**2 + ( (X1-X2)*6.5/40 )**2 + ( (Y1-Y2)*6.5/40 )**2 ) )
			data[0].append( np.sqrt( ( (Z1-Z2)*2.5 )**2 ) )
			data[1].append( trow.times )
		
	ax1.plot( data[1], data[0], 'ob', ms = 3, mew = 0 )

	data = np.transpose( np.array(data) )
	# for i,val in enumerate(average[0]):
	# 	average[1][i] += np.mean([i[0] for i in data if i[1] < val and i[1]>=(val-6)])
	oldaverage += np.nanmean([i[0] for i in data if i[1] >20])
# average = np.array(average)
# average[1] /= 8
# print(average)
	# print(data)
print(oldaverage/8)





### plot XYdiameter

worms = [ i[-10:-7] for i in worms ]

# average = [[3,9,15,21,27,33],[0,0,0,0,0,0]]
oldaverage = 0
for worm in worms:
	data = np.loadtxt( open( path + '\\XYdiameter_' + worm  + '.txt', 'rb' ) )
	times = data[1:,2] / 60.
	lengths = data[1:,3] * 6.5 / data[1:,4]
	print(times, lengths)
	ax1.plot( times, lengths, 'or',
		ms = 3,
		fillstyle='full',
		mew = 0 )

	data = np.transpose( [lengths,times] )
	# for i,val in enumerate(average[0]):
	# 	average[1][i] += np.nanmean([i[0] for i in data if i[1] < val and i[1]>=(val-6)])
	oldaverage += np.nanmean([i[0] for i in data if i[1] >20])

# average = np.array(average)
# average[1] /= 8
# print(average)
print(oldaverage/8)

plt.show()
