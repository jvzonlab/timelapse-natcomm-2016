import glob
from timelapseFun import *
from tifffile import *
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.stats import gaussian_kde

import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def straighten( path, channels = ['LED','488nm'], worms = 'all' ):

	if worms == 'all':
		folders = [ i for i in os.listdir(path) if os.path.isdir(path+i) ]
		foldersToKeep = [ ( i.endswith('downsized') ) for i in folders ]
		worms = [ i for idx, i in enumerate(folders) if foldersToKeep[idx] ]
		worms.sort()

	elif isinstance(worms,list):
		worms = [ i+'_downsized' for i in worms ]
		# worms = [ i for i in worms ]

	for w in worms:
		print( path, w )
		# path = 'Y:\\Images\\150821_lin42_mlt10_250x250x20\\'+w+'_downsized'
		wormpath = path + w

		# load worm data
		loadPath = path + 'worm' + w.split('_')[0] + '.pickle'
		df = pickle.load(open(loadPath,'rb'))
		bodyData = df.ix[ df.rowtype == 'body', ['outline', 'times'] ]

		times = bodyData.times.values#ix[ (bodyData.times>=0) & (len(bodyData.outline)>2), 'times' ]

		for channel in channels:
			flist = glob.glob(wormpath+'\\*' + channel + '.tif')
			flist.sort()
			movie = []
			_min = 2**16
			_max = 0

			for f in flist[:]:
				print(f)
				imgs = loadstack( f )
				if channel == 'LED':
					img = np.mean(imgs,0).astype(np.uint16)
				else:
					img = np.max(imgs,0).astype(np.uint16)
				movie.append( img )

				_min = np.min(img) * (np.min(img)<_min) + _min * (np.min(img)>=_min)
				_max = np.max(img) * (np.max(img)>_max) + _max * (np.max(img)<=_max)

			### RESIZE ALL THE IMAGES BASED ON THE MAXIMUM LENGTH
			movieFinal = np.zeros( ( len(movie), 2 * movie[0].shape[0], 2 * movie[0].shape[1] ) ).astype( np.uint8 )
			movieFinalWithTime = np.zeros( ( len(movie), 2 * movie[0].shape[0], 2 * movie[0].shape[1] ) ).astype( np.uint8 )

			for idx, img in enumerate( movie ):
				print(idx)
				img = ( 2**8 - 1. ) * ( img - _min ) / ( _max - _min )

				f = ip.interp2d(np.arange(img.shape[0]), np.arange(img.shape[1]), img, kind='linear')
				imgnew = f(np.linspace(0, img.shape[0],2*img.shape[0]),np.linspace(0, img.shape[0],2*img.shape[0]))

				movieFinal[ idx ] = imgnew.astype( np.uint8 )

				imgpil = Image.fromarray( movieFinal[idx], 'L' )

				font = ImageFont.truetype( "calibri.ttf", 60 )
				draw = ImageDraw.Draw(imgpil)
				draw.text((0,0),'%d h' % np.floor(times[idx]),fill='white',font=font)

				movieFinalWithTime[idx] = np.asarray( imgpil )


			imsave( wormpath + '\\' + channel + 'movie.tif', np.array(movieFinal) )
			# imsave( wormpath + '\\' + channel + 'movieWithTime.tif', np.array( movieFinalWithTime ) )


if __name__ == '__main__':

	path = 'X:\\Images\\150914_mlt10_250x250x20\\'

	straighten( path, channels = ['LED', '488nm'], worms = [ 'C01' ] )
