# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:52:19 2015

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
from scipy.stats import gaussian_kde

import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

################### setup the figure ###################################
# fig, ax = plt.subplots(figsize=(1,1*15.18/20.32))
fig, ax = plt.subplots(figsize=(2,2*35.677/42.84))
# fig.subplots_adjust(left=0.15, right=.88, top=.97, bottom=0.15)
fig.subplots_adjust(left=0., right=1., top=1., bottom=0.)

#create axes and plot
ax2 = ax.twinx()
ax2.set_ylim(0,1.1)
ax2.set_xlim(0,50)
ax2.set_yticks(np.arange(0.2,1.01,.2))

# set the right axes to red
for tl in ax2.get_yticklabels():
   tl.set_color('b')
   tl.set_fontsize(15)

# set the left axes to be black and set ylim
for tl in ax.get_yticklabels():
   tl.set_color('k')
   tl.set_fontsize(15)
for tl in ax.get_xticklabels():
   tl.set_fontsize(15)
ax.set_xlim(0,50)
ax.set_ylim(0,1.2)
ax.set_yticks(np.arange(0.2,1.01,.2))
ax.set_xticks(np.arange(10,50,10))

ax.set_xlabel('Hours after hatching',fontsize=15)
ax.set_ylabel('Body length [mm]',fontsize=15)
ax2.text(55,.6,'Fraction',fontsize = 15,rotation=-90,color = 'b')

# horixontal line for 100%
ax2.plot( [0, 50], [1, 1], 'k--' )

correctForStrain = False
correctForSingleWorms = False

############################### 290x290x25

path = 'Y:\\Images\\150810_N2_290x290x25\\'

flist = glob.glob(path+'wormC*.pickle')
flist.sort()

color = 'red'
ax.text( 2, 1., '290x290x25, N=%d'%len(flist), color = color, fontsize = 15)

### plot all the length

allLength = []
for w in flist:
    print(w)
    df = pickle.load( open( w, 'rb' ) )

    if correctForStrain:
      df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
    if correctForSingleWorms:
      df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

    mag = df.ix[df.rowtype == 'param', 'magnification'].values[0]
    com = df.ix[df.rowtype == 'param', 'compression'].values[0]
#   print(mag,com)

    body = df.ix[(df.rowtype=='body') & (df.tidx>=0), ['times', 'tidx', 'spline'] ]

    length = []
    for spline in body.spline:
        if len(spline) > 1:
            length.append( len(spline) * com * 6.5 / ( mag * 1000 ) )
        else:
            length.append(np.nan)
    body['length'] = length

    bins = np.arange(0,50,.5)
    groups = body.groupby( pd.cut( body.times, bins ) )
    allLength.append( groups.mean().length.values )

    ax.plot( bins[:-1], groups.mean().length.values, '-', color = color, alpha = .1 )

allLength = np.array(allLength)
ax.plot( bins[:-1], np.nanmean(allLength, 0), color = color, lw = 2 )

### plot the ecdysis 

ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )
totLeth = []

for w in flist:
    worm = w[-10:-7]
    
    df = pickle.load( open( w, 'rb' ) )

    if correctForStrain:
      df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
    if correctForSingleWorms:
      df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

    # load ecdysis data
    index = np.where( ecd[:,0] == float(worm[1:]) )
    mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
    lethtidx = ecd[ index, 2:6 ][0][0] - mintp - 1

    lethTp = [np.nan,np.nan,np.nan,np.nan]
    for idx, tidx in enumerate( lethtidx ):
      if tidx > 0:
        lethTp[idx] = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]

    totLeth.append(lethTp)

totLeth = np.array(totLeth)    

for idx in np.arange(4):
    times = totLeth[:,idx]
    times = times[~np.isnan(times)]
    # times = ( ecd[:,idx+2] - ecd[:,1] ) /6.
    # times = np.array( [ i for i in times if i>0 ] )
    
    _min = np.floor(np.min(times)-5)
    _max = np.floor(np.max(times)+5)

    # distr, bins = np.histogram( times, bins = (_max - _min)*2, range = [_min,_max] )
    # ax2.fill( bins[:-1], distr / np.sum(distr), ec=None, fc=color, alpha=0.4)

    weights = np.ones_like(times)/len(times)
    n, bins, rectangles = ax2.hist( times, 100, range = (0,50), weights=weights, color = color, alpha = .6, rwidth = .8, lw=0)

    # kernel = gaussian_kde( times )
    # ax2.fill( bins, kernel(bins)*(bins[1]-bins[0]), ec=None, fc=color, alpha=0.4)    



############################### 250x250x20

path = 'Y:\\Images\\150731_N2_250x250x20\\'

flist = glob.glob(path+'*.pickle')
flist.sort()

color = 'blue'
# ax.text( 2, .92, '250x250x20, N=%d'%len(flist), color = color, fontsize = 15)
ax.text( 2, 1.02, 'N=%d'%len(flist), color = 'black', fontsize = 15)

### plot all the length

allLength = []
for w in flist:
   print(w)
   df = pickle.load( open( w, 'rb' ) )

   if correctForStrain:
     df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
   if correctForSingleWorms:
     df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

   mag = df.ix[df.rowtype == 'param', 'magnification'].values[0]
   com = df.ix[df.rowtype == 'param', 'compression'].values[0]

   body = df.ix[(df.rowtype=='body') & (df.times>=0), ['times', 'tidx', 'spline'] ]

   length = []
   for spline in body.spline:
       if len(spline) > 1:
           length.append( len(spline) * com * 6.5 / ( mag * 1000 ) )
       else:
           length.append(np.nan)
   body['length'] = length

   body = body.ix[ pd.notnull(body.length) ]

   bins = np.arange(0,50,.5)
   groups = body.groupby( pd.cut( body.times, bins ) )
   allLength.append( groups.mean().length.values )

   ax.plot( bins[:-1], groups.mean().length.values, '-', color = color, alpha = .2 )

allLength = np.array(allLength)
ax.plot( bins[:-1], np.nanmean(allLength, 0), color = color, lw = 2 )

### highlight images

selectedWorm = path + 'wormC12.pickle'
selectedTimes = [8,16,24,32,40]

df = pickle.load( open( selectedWorm, 'rb' ) )

if correctForStrain:
 df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
if correctForSingleWorms:
 df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

mag = df.ix[df.rowtype == 'param', 'magnification'].values[0]
com = df.ix[df.rowtype == 'param', 'compression'].values[0]

body = df.ix[(df.rowtype=='body') & (df.times>=0), ['times', 'tidx', 'spline'] ]

length = []
for spline in body.spline:
   if len(spline) > 1:
       length.append( len(spline) * com * 6.5 / ( mag * 1000 ) )
   else:
       length.append(np.nan)
body['length'] = length

body = body.ix[ pd.notnull(body.length) ]

selectedPoints = []
for t in selectedTimes:
  selectedPoints.append([body.ix[ ( (body.times - t) < .5 ) & ( (body.times - t) > 0 ), [ 'times' ] ].values[0], 
                        body.ix[ ( (body.times - t) < .5 ) & ( (body.times - t) > 0 ), [ 'length' ] ].values[0]]) 
# print(selectedPoints)
selectedPoints = np.array(selectedPoints)
ax.plot( selectedPoints[:,0], selectedPoints[:,1], 'or',
         fillstyle='full',
         mec='red' )


### plot the ecdysis 

ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )
totLeth = []

for w in flist:
    worm = w[-10:-7]
    
    df = pickle.load( open( w, 'rb' ) )

    if correctForStrain:
      df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
    if correctForSingleWorms:
      df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

    # load ecdysis data
    index = np.where( ecd[:,0] == float(worm[1:]) )
    mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
    lethtidx = ecd[ index, 2:6 ][0][0] - mintp - 1

    lethTp = [np.nan,np.nan,np.nan,np.nan]
    for idx, tidx in enumerate( lethtidx ):
      if tidx > 0:
        lethTp[idx] = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]

    totLeth.append(lethTp)

totLeth = np.array(totLeth)    

for idx in np.arange(4):
    times = totLeth[:,idx]
    times = times[~np.isnan(times)]
    # times = ( ecd[:,idx+2] - ecd[:,1] ) /6.
    # times = np.array( [ i for i in times if i>0 ] )
    
    _min = np.floor(np.min(times)-5)
    _max = np.floor(np.max(times)+5)

    # distr, bins = np.histogram( times, bins = (_max - _min)*2, range = [_min,_max] )
    # ax2.fill( bins[:-1], distr / np.sum(distr), ec=None, fc=color, alpha=0.4)

    weights = np.ones_like(times)/len(times)
    n, bins, rectangles = ax.hist( times, 100, range = (0,50), weights=weights, color = color, alpha = .6, rwidth = .8, lw=0)

    # kernel = gaussian_kde( times )
    # ax2.fill( bins, kernel(bins)*(bins[1]-bins[0]), ec=None, fc=color, alpha=0.4)    

# ############################# N2 ON PLATES - small OP50 area

# path = 'Y:\\Images\\160513_N2_on_plates'
# worms = ['w1_1','w1_2','w1_3','w1_4','w2_4','w3_1','w3_2','w4_1','w4_2','w4_3','w4_4']

# for worm in worms:
#   data = np.loadtxt( open( path + '\\' + worm  + '.txt', 'rb' ) )
#   times = data[1:,2] / 60.
#   lengths = data[1:,3] * 6.5 / data[1:,4] / 1000.
#   print(times, lengths)
#   ax.plot( times, lengths, 'ow',
#            fillstyle='full',
#            mec='black' )

# ############################# N2 ON PLATES - OP50 on all plate

# path = 'Y:\\Images\\160520_N2_on_plates\\20deg_Bseries'
# worms = ['wB01','wB02','wB03','wB04','wB05','wB06','wB07','wB08','wB09','wB10','wB11','wB12','wB13','wB14','wB15']

# for worm in worms:
#   data = np.loadtxt( open( path + '\\' + worm  + '.txt', 'rb' ) )
#   times = data[1:,2] / 60.
#   lengths = data[1:,3] * 6.5 / data[1:,4] / 1000.
#   print(times, lengths)
#   ax.plot( times, lengths, 'ob',
#            fillstyle='full',
#            mec='black' )

############################# N2 ON PLATES - OP50 on all plate at 22 degrees

path = 'Y:\\Images\\160520_N2_on_plates\\22deg_Aseries'
worms = ['wA01','wA02','wA04','wA05','wA06','wA07','wA08','wA09','wA12','wA15']

for worm in worms:
  data = np.loadtxt( open( path + '\\' + worm  + '.txt', 'rb' ) )
  times = data[1:,2] / 60.
  lengths = data[1:,3] * 6.5 / data[1:,4] / 1000.
  print(times, lengths)
  ax.plot( times, lengths, 'ow',
           fillstyle='full',
           mec='black' )

# ############################### mlt10

# path = 'Y:\\Images\\150914_mlt10_250x250x20\\'

# flist = glob.glob(path+'wormC*.pickle')
# flist.sort()

# color = 'yellow'
# ax.text( 2, 1., '290x290x25, N=%d'%len(flist), color = color, fontsize = 15)

# ### plot all the length

# allLength = []
# for w in flist:
#     print(w)
#     df = pickle.load( open( w, 'rb' ) )

#     if correctForStrain:
#       df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
#     if correctForSingleWorms:
#       df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

#     mag = 40
#     com = 8
# #   print(mag,com)

#     body = df.ix[(df.rowtype=='body') & (df.tidx>=0), ['times', 'tidx', 'spline'] ]

#     length = []
#     for spline in body.spline:
#         if len(spline) > 1:
#             length.append( len(spline) * com * 6.5 / ( mag * 1000 ) )
#         else:
#             length.append(np.nan)
#     body['length'] = length

#     bins = np.arange(0,50,.5)
#     groups = body.groupby( pd.cut( body.times, bins ) )
#     allLength.append( groups.mean().length.values )

#     ax.plot( bins[:-1], groups.mean().length.values, '-', color = color, alpha = .1 )

# allLength = np.array(allLength)
# ax.plot( bins[:-1], np.nanmean(allLength, 0), color = color, lw = 2 )

# ### plot the ecdysis 

# ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )
# totLeth = []

# for w in flist:
#     worm = w[-10:-7]
    
#     df = pickle.load( open( w, 'rb' ) )

#     if correctForStrain:
#       df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
#     if correctForSingleWorms:
#       df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

#     # load ecdysis data
#     index = np.where( ecd[:,0] == float(worm[1:]) )
#     mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
#     lethtidx = ecd[ index, 2:6 ][0][0] - mintp - 1

#     lethTp = [np.nan,np.nan,np.nan,np.nan]
#     for idx, tidx in enumerate( lethtidx ):
#       if tidx > 0:
#         lethTp[idx] = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]

#     totLeth.append(lethTp)

# totLeth = np.array(totLeth)    

# for idx in np.arange(4):
#     times = totLeth[:,idx]
#     times = times[~np.isnan(times)]
    
#     _min = np.floor(np.min(times)-5)
#     _max = np.floor(np.max(times)+5)
    
#     distr, bins = np.histogram( times, bins = (_max - _min)*1.5, range = [_min,_max] )
    
#     kernel = gaussian_kde( times )
    
#     ax2.fill( bins, kernel(bins)*(bins[1]-bins[0]), ec=None, fc=color, alpha=0.4)    


# ############################### wrt2

# path = 'Z:\\150824_JVZ20_wrt2GFP_250x250x20\\'

# flist = glob.glob(path+'wormC*.pickle')
# flist.sort()

# color = 'green'
# ax.text( 2, 1., '290x290x25, N=%d'%len(flist), color = color, fontsize = 15)

# ### plot all the length

# allLength = []
# for w in flist:
#     print(w)
#     df = pickle.load( open( w, 'rb' ) )

#     if correctForStrain:
#       df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
#     if correctForSingleWorms:
#       df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

#     mag = 40
#     com = 8
# #   print(mag,com)

#     body = df.ix[(df.rowtype=='body') & (df.tidx>=0), ['times', 'tidx', 'spline'] ]

#     length = []
#     for spline in body.spline:
#         if len(spline) > 1:
#             length.append( len(spline) * com * 6.5 / ( mag * 1000 ) )
#         else:
#             length.append(np.nan)
#     body['length'] = length

#     bins = np.arange(0,50,.5)
#     groups = body.groupby( pd.cut( body.times, bins ) )
#     allLength.append( groups.mean().length.values )

#     ax.plot( bins[:-1], groups.mean().length.values, '-', color = color, alpha = .1 )

# allLength = np.array(allLength)
# ax.plot( bins[:-1], np.nanmean(allLength, 0), color = color, lw = 2 )

# ### plot the ecdysis 

# ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )
# totLeth = []

# for w in flist:
#     worm = w[-10:-7]
    
#     df = pickle.load( open( w, 'rb' ) )

#     if correctForStrain:
#       df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
#     if correctForSingleWorms:
#       df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

#     # load ecdysis data
#     index = np.where( ecd[:,0] == float(worm[1:]) )
#     mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
#     lethtidx = ecd[ index, 2:6 ][0][0] - mintp - 1

#     lethTp = [np.nan,np.nan,np.nan,np.nan]
#     for idx, tidx in enumerate( lethtidx ):
#       if tidx > 0:
#         lethTp[idx] = df.ix[(df.rowtype=='body')&(df.tidx==tidx),'times'].values[0]

#     totLeth.append(lethTp)

# totLeth = np.array(totLeth)    

# for idx in np.arange(4):
#     times = totLeth[:,idx]
#     times = times[~np.isnan(times)]
    
#     _min = np.floor(np.min(times)-5)
#     _max = np.floor(np.max(times)+5)
    
#     distr, bins = np.histogram( times, bins = (_max - _min)*1.5, range = [_min,_max] )
    
#     kernel = gaussian_kde( times )
    
#     ax2.fill( bins, kernel(bins)*(bins[1]-bins[0]), ec=None, fc=color, alpha=0.4)    







plt.show()