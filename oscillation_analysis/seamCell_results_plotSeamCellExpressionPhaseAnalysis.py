# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:30:41 2015

@author: gritti
"""

'''

TO BE IMPLEMENTED!!!

'''
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
from skimage import morphology
from scipy.signal import detrend
mpl.rcParams['pdf.fonttype'] = 42
#from uiget_dir import *
#import psutil

def loadData( path, worm, cell1, cell2, cside ):

    df = pickle.load(open(path+'worm'+worm+'.pickle','rb'))

    cellmask = df.rowtype=='cell'
    cellsData = df.ix[ cellmask ]

    cell1Data = cellsData.ix[ ( cellsData.cname.str.startswith( cell1 ) ) & ( cellsData.cside == cside ), ['cname','wrt2exp','tidx','times'] ].groupby(cellsData.tidx).mean()
    cell1Data.cname = cell1
    cell2Data = cellsData.ix[ ( cellsData.cname.str.startswith( cell2 ) ) & ( cellsData.cside == cside ), ['cname','wrt2exp','tidx','times'] ].groupby(cellsData.tidx).mean()
    cell2Data.cname = cell2

    cell1Data.exp = detrend( cell1Data.wrt2exp )
    cell2Data.exp = detrend( cell2Data.wrt2exp )

    dT = 20
    # filter data for cell 1
    times = np.linspace( np.min(cell1Data.times), np.max(cell1Data.times), 60 / dT * ( np.max(cell1Data.times) - np.min(cell1Data.times) ) )
    interpolation = ip.interp1d(cell1Data.times,cell1Data.wrt2exp, kind='linear')
    interpexp = interpolation(times)

    c1Data = pd.DataFrame( {
        'cname': cell1,
        'tidx': np.arange(len(interpexp)),
        'times': times,
        'exp': detrend( filters.gaussian_filter1d( interpexp, 20/dT ) )
        } )

    # filter data for cell 1
    times = np.linspace( np.min(cell2Data.times), np.max(cell2Data.times), 60 / dT * ( np.max(cell2Data.times) - np.min(cell2Data.times) ) )
    interpolation = ip.interp1d(cell2Data.times,cell2Data.wrt2exp, kind='linear')
    interpexp = interpolation(times)

    c2Data = pd.DataFrame( {
        'cname': cell2,
        'tidx': np.arange(len(interpexp)),
        'times': times,
        'exp': detrend( filters.gaussian_filter1d( interpexp, 20/dT ) )
        } )

    return ( cell1Data, cell2Data, c1Data, c2Data )

def getDeviation( c1, c2 ):

    # std of distance from diagonal
    d = []
    for c in zip(c1.exp,c2.exp):
        p = np.array([c[0],c[1]])
        alpha = np.pi/4. - np.arctan2(c[1],c[0])
        # print(p,alpha)

        d.append(np.abs(np.sqrt(np.sum(p*p))*np.sin(alpha)))

    # x = c1.exp
    # y = c2.exp
    # d = np.mean( (x-np.mean(x))*(y-np.mean(y)) ) / (np.std(x)*np.std(y))
    return np.mean(d)





path = 'X:\\Images\\150824_JVZ20_wrt2GFP_250x250x20\\'
worm = 'C10'

cells1 = ['b.','6.','c.','5.','4.','1.','2.']
cell2 = '3.'
cside = 'L'

fig = plt.figure(figsize = (5.8,5.8))
# plt.title(cell1+'-'+cell2)
ax = fig.add_subplot(111)
ax.set_xlim((-10000,10000))
ax.set_ylim((-10000,10000))

colors = np.array( [cm.brg(int(i)) for i in np.arange(len(cells1))*254./(len(cells1)-1)] )
colors = colors[::-1]

### find deviation from V3

dev = []
for idx, cell1 in enumerate( cells1 ):
    c1Raw, c2Raw, c1Filt, c2Filt = loadData(path, worm, cell1, cell2, cside)
    dev.append( getDeviation( c1Raw, c2Raw ) )

### sort the data according to deviation from V3

cells1sorted = [x for (y,x) in sorted(zip(dev,cells1))][::-1]
devsorted = [y for (y,x) in sorted(zip(dev,cells1))][::-1]

for idx, cell1 in enumerate(cells1sorted):
    print(cell1,devsorted[idx])
    c1Raw, c2Raw, c1Filt, c2Filt = loadData(path, worm, cell1, cell2, cside)
    ax.plot( c1Raw.exp, c2Raw.exp, 'o', label = cell1 + '-' + str(int(devsorted[idx])), color = colors[idx] )

plt.legend(loc=0)
plt.show()
