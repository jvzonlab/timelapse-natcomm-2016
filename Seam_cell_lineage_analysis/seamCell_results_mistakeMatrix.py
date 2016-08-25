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


data = np.zeros((4,10))

### TOTAL
# data[1,3] = 9/48
# data[1,4] = 1/48
# data[1,8] = 1/48
# data[2,1] = 6/16
# data[2,2] = 2/16
# data[2,3] = 1/32
# data[2,6] = 2/32
# data[3,1] = 7/16
# data[3,2] = 2/16
# data[3,3] = 3/32
# data[3,4] = 7/32
# data[3,5] = 5/32
# data[3,6] = 5/32
# data[3,7] = 3/16
# data[3,8] = 17/32

### AS->SYM
data[1,3] = 4/48
data[1,4] = 1/48
data[1,8] = 1/48
data[2,1] = 6/16
data[2,2] = 2/16
data[2,3] = 1/32
data[2,6] = 0/32
data[3,1] = 5/16
data[3,2] = 2/16
data[3,3] = 3/32
data[3,4] = 7/32
data[3,5] = 5/32
data[3,6] = 5/32
data[3,7] = 3/16
data[3,8] = 17/32

# ### NO DIV
data[1,3] = 5/48
data[1,4] = 0/48
data[1,8] = 0/48
data[2,1] = 0/16
data[2,2] = 0/16
data[2,3] = 0/32
data[2,6] = 2/32
data[3,1] = 2/16
data[3,2] = 0/16
data[3,3] = 0/32
data[3,4] = 0/32
data[3,5] = 0/32
data[3,6] = 0/32
data[3,7] = 0/16
data[3,8] = 0/32




nullfmt = NullFormatter()

fig = plt.figure(figsize = (5.8,5.8))

left, width = 0.1, 0.6
bottom, height = 0.3, 0.24
bottom_h = bottom  - width/10 - 0.01
left_h = left + width + .01
left_hh = left#left_h + height/4 + 0.01
bottom_hh = bottom_h - width/10 - 0.01

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, width/10]
rect_histy = [left_h, bottom, height/4, height]
rect_cmap = [left_hh, bottom_hh, width, width/10]

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
axCmap = plt.axes(rect_cmap)

axScatter.xaxis.tick_top()
axCmap.yaxis.tick_right()

axScatter.xaxis.set_ticks_position('none')
axScatter.yaxis.set_ticks_position('none')
axHistx.xaxis.set_ticks_position('none')
axHistx.yaxis.set_ticks_position('none')
axHisty.xaxis.set_ticks_position('none')
axHisty.yaxis.set_ticks_position('none')
axCmap.xaxis.set_ticks_position('none')


# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHistx.yaxis.set_major_formatter(nullfmt)
axHisty.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)
axCmap.yaxis.set_major_formatter(nullfmt)

# invert Y hist axes
axHisty.invert_yaxis()

maxColor = .5

cax = axScatter.imshow(data,interpolation='nearest', cmap=cm.Blues, clim = [0,maxColor])

axCmap.set_xticks([0,128,255])
axCmap.set_xticklabels([0,maxColor/2,maxColor],fontsize=15)

# axHistx.set_clim(0,1.)
# axHisty.set_clim(0,1.)

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))[::-1]
cmapimg = axCmap.imshow(gradient, aspect='auto', cmap=cm.Blues)

axHistx.imshow(np.vstack((np.mean(data,0),np.mean(data,0))), cmap=cm.Blues,interpolation='nearest', aspect='auto',clim = [0,maxColor])

axHisty.imshow(np.vstack((np.mean(data,1),np.mean(data,1))).T, cmap=cm.Blues,interpolation='nearest', aspect='auto',clim = [0,maxColor])


axHistx.set_xticks(np.linspace(0,8,9))
minorLocator = AutoMinorLocator(2)
axHistx.xaxis.set_minor_locator(minorLocator)
axHistx.grid(which='minor')

axHisty.set_yticks(np.linspace(0,3,4))
minorLocator = AutoMinorLocator(2)
axHisty.yaxis.set_minor_locator(minorLocator)
axHisty.grid(which='minor')

axScatter.set_xticks([0,1,2,3,4,5,6,7,8,9])
axScatter.set_yticks([0,1,2,3])
axScatter.set_xticklabels(['H0','H1','H2','V1','V2','V3','V4','V5','V6','T'], minor=False, fontsize = 15)
axScatter.set_yticklabels(['L1','L2','L3','L4'], fontsize=15)
minorLocator = AutoMinorLocator(2)
axScatter.yaxis.set_minor_locator(minorLocator)
axScatter.xaxis.set_minor_locator(minorLocator)
axScatter.grid(which='minor')

plt.show()
