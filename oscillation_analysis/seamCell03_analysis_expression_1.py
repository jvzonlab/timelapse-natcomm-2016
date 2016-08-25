# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:22:21 2015

@author: gritti
"""


import glob
from timelapseFun import *
from matplotlib import cm
from skimage import morphology, filter
import numpy as np
import pickle


# path = 'Y:\\Images\\150617_JVZ20_wrt2GFP\\'
path = 'X:\\Images\\150824_JVZ20_wrt2GFP_250x250x20\\'
worm = 'C22'

color = np.array( [cm.gist_rainbow(int(i))[:3] for i in np.arange(4)*200./3.] )

df = pickle.load( open(path + 'worm' + worm + '.pickle','rb') )

flist = glob.glob(path + worm + '_straighten\\straight*_488nm.tif')
flist.sort()
flistnew = ['' for i in df.ix[(df.rowtype=='body')&(df.tidx>=0),'tidx']]

firsttidx = int(flist[0].split('_')[-2][-3:]) - 1
for f in flist:
    tidx = int(f.split('_')[-2][-3:])-1
    flistnew[tidx-firsttidx] = f

for i in np.arange(np.min(df.ix[df.rowtype=='cell','tidx'])):
    flistnew.insert(0,'')
    
cellrowmask = df.rowtype == 'cell'
bodyrowmask = df.rowtype == 'body'

radius = 10

for tidx, f in enumerate( flistnew ):
    print(f)
    
    tidxmask = df.tidx == tidx
    cells = df.ix[ tidxmask & cellrowmask ]
    
    if len(cells) > 0:
        
        imgs = loadstack( f )
        
        for idx, cell in cells.iterrows():
            
            cimg = imgs[ cell.cZpos, 
                        np.clip(cell.cYpos-radius,0,None):cell.cYpos+radius+1, 
                        np.clip(cell.cXpos-radius,0,None):cell.cXpos+radius+1 ]

            cmask = cimg > filter.threshold_otsu(cimg)
            signal = np.sum( cimg * cmask ) / np.sum( cmask )
#            fig = plt.figure()
#            axarr = fig.add_subplot(121)
#            axarr.imshow(cimg,cmap='gray',interpolation='nearest')
#            axarr = fig.add_subplot(122)
#            axarr.imshow(cmask,interpolation='nearest')

            df.ix[idx,'wrt2exp'] = signal

pickle.dump( df, open( path + 'worm'+worm+'.pickle', 'wb' ), protocol = 2 )
             

