import glob
from timelapseFun import *
from tifffile import *
import numpy as np
from PIL import Image
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.stats import gaussian_kde

import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42



def drawKymographOneWorm( path='Y:\\Images\\150617_JVZ20_wrt2GFP\\', worm='C01' ):

	# setup the figure
	fig = plt.figure(figsize = (5,5))
	ax = fig.add_subplot(111)
	ax.axis('off')
	fig.subplots_adjust(left=0., right=1., top=1., bottom=0.)

	# LOAD DATAFRAME and 
	df = pickle.load( open( path+'worm'+worm+'.pickle', 'rb' ) )
	print(df.keys())
	paramdata = df.ix[ df.rowtype == 'param' ]
	bodydata = df.ix[ df.rowtype == 'body', ['tidx','times','spline','length'] ]
	cellsdata = df.ix[ df.rowtype == 'cell' ]
	print(2048/paramdata.compression.values[0])

	# compute image lengths, calculate max length and the center of the kymograph
	bodydata.length = [ (len(i)>1) * (len(interp_spline(2048/paramdata.compression.values[0]*i))-2)+(len(i)==1)*0 for i in bodydata.spline.values]
	bodydata = bodydata.ix[ bodydata.length > 0 ].reset_index(True)
	maxLength = np.max( bodydata.length )
	center = int( ( maxLength + 1 )/2 )
	print(maxLength)

	# create blank kymograph
	kymogr = np.zeros( ( len(bodydata.length), maxLength ) ).astype( np.float )

	# read in the file list of images
	wormpath = path + worm + '_straighten'
	flist = glob.glob(wormpath+'\\*488nm.tif')
	flist.sort()

	# start creating the kymograph
	_min = 2**16
	_max = 0
	for idx, f in enumerate( flist[:] ):
		# load the image and perform mean in uint16
		print(f)
		imgs = loadstack( f )
		img = np.mean( np.max(imgs,0), 0).astype('uint16')

		# update the max and min value of the kymograph
		_min = np.min(img) * (np.min(img)<_min) + _min * (np.min(img)>=_min)
		_max = np.max(img) * (np.max(img)>_max) + _max * (np.max(img)<=_max)

		# place the new line in the kymograph
		kymogr[ idx : ( idx + 1 ), int( center - img.shape[0]/2 ) : int( center + img.shape[0]/2 ) ] = img.astype( np.uint16 )

	# adjust brightness and contrast
	kymogr = (2**16-1)*(kymogr - np.min(kymogr)) / (np.max(kymogr)-np.min(kymogr)) 
	ax.imshow( kymogr, cmap='gray', interpolation = 'nearest', aspect = 'auto' )

	# plot worm head and tail
	ax.plot( (maxLength - bodydata.length)/2, np.arange(len(bodydata.length)), '--w', lw=1, dashes = [2,2] )
	ax.plot( maxLength - (maxLength - bodydata.length)/2, np.arange(len(bodydata.length)), '--w', lw = 1, dashes = [2,2] )

	## PLOT ECDYSIS

	# load txt file with ecdysis indexes
	ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )

	# load ecdysis indexes for the worm
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - mintp - 1

	# find ecdysis times
	lethTp = [np.nan,np.nan,np.nan,np.nan]
	for idx, tidx in enumerate( lethtidx ):
	  if tidx > 0:
	    lethTp[idx] = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]
	print('lethargus times:', lethTp)

	# find the row in the kymograph and plot dash line
	for leth in lethTp:
		row = bodydata.ix[ bodydata.times == leth ].index.values[0] + 1
		ax.plot( [0,maxLength], [row,row], '--w', dashes = [2,2], lw = 1 )
		print(row)

	## draw color bars at 25,50 and 75%
	lines = [[],[],[]]
	for bl in bodydata.length:
		lines[0].append(bl*25/100+(maxLength-bl)/2)
		lines[1].append(bl*50/100+(maxLength-bl)/2)
		lines[2].append(bl*75/100+(maxLength-bl)/2)
	ax.plot( lines[0], np.arange(len(bodydata.length)), '-y', lw=2, alpha = .7 )
	ax.plot( lines[1], np.arange(len(bodydata.length)), '-r', lw=2, alpha = .7 )
	ax.plot( lines[2], np.arange(len(bodydata.length)), '-c', lw=2, alpha = .7 )


	ax.set_xlim([0,maxLength])
	ax.set_ylim([0,len(bodydata.length)])
	ax.invert_yaxis()
	plt.show()

if __name__ == '__main__':

	path = 'Y:\\Images\\150914_mlt10_250x250x20\\'
	worm = 'C10'

	drawKymographOneWorm( path, worm,  )
