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
import gc
mpl.rcParams['pdf.fonttype'] = 42

def plotExpressionOneWorm( axes, gene, path='Y:\\Images\\150617_JVZ20_wrt2GFP\\', worm='C01', correctForStrain = False, correctForSingleWorms = False,
	plotDots = True, plotMean = True, plotSmooth = True, plotLeth = True, plotNSeam = True, color = [0,0,0,1], colorDots = [0.8,0.8,0.8], alpha = .5, 
	lineage = 'all', lw = 3, shift = None ):

	# load the data    
	loadPath = path + 'worm' + worm + '.pickle'
	df = pickle.load(open(loadPath,'rb'))

	if correctForStrain:
		df.times *= df.ix[ df.rowtype == 'param', 'strainScaleTimeFactor' ].values[0]
	if correctForSingleWorms:
		df.times *= df.ix[ df.rowtype == 'param', 'scaleTimeFactor' ].values[0]

	if gene=='mlt10':
		df.ix[:,'exp'] = df.exp / 1000
		exp = df.ix[df.rowtype == 'body', ['tidx','times','exp']]
	elif gene=='wrt2':
		df.ix[:,'exp'] = df.wrt2exp / 1000
		exp = df.ix[ df.rowtype == 'cell', ['tidx','times','cname','cside','exp'] ]
		vfilter = exp.cname.str.startswith('1.') | exp.cname.str.startswith('2.') | exp.cname.str.startswith('3.') | exp.cname.str.startswith('4.') | exp.cname.str.startswith('5.')
		exp = exp.ix[ vfilter ]
		if lineage != 'all':
			vfilter = exp.cname.str.startswith(lineage)
			exp = exp.ix[ vfilter ]
	# print(exp)

	# plot the gray dots for single cells
	if plotDots:
		if lineage == 'all':
		    axes.plot(exp.times,exp.exp,'o',color = colorDots, alpha = alpha, mec = None, mew = 0.0, markersize = 5.)
		else:
			axes.plot( exp.times, exp.exp + shift, 'o', color = color, mec = None, mew = 0.0, markersize = 3. )

	## PLOT THE MEAN
	# plot the mean expression level binned every 1/N hour
	bins = np.arange(0,50,.5)
	groups = exp.groupby(pd.cut(exp.times,bins))

	if plotMean:
	    pl1, = axes.plot(groups.mean().times[pd.notnull(groups.mean().exp)],
	        groups.mean().exp[pd.notnull(groups.mean().exp)],'-',lw=2,color = color)

	## PLOT THE FILTERED DATA

	if plotSmooth:
	    # group data by time and remove missing datapoint
	    groups = exp.groupby(exp.times)
	    data = groups.mean().ix[ (groups.mean().tidx >= 0) ]
	    data = data.ix[pd.notnull(data.exp)]
	    data.ix[:,'times'] = data.index
	    data = data.reset_index(drop=True)

	    # compute finer 1d interpolation
	    times = np.linspace( np.min(data.times), np.max(data.times), 60 * ( np.max(data.times) - np.min(data.times) ) )
	    interpolation = ip.interp1d(data.times,data.exp, kind='linear')
	    interpexp = interpolation(times)

	    # compute gaussian filter of data
	    filtdata = pd.DataFrame( {
	        'times': times,
	        'filtexp': filters.gaussian_filter1d(interpexp,60)
	        } )
	    pl1, = axes.plot(filtdata.times,
	            filtdata.filtexp,'-',lw=lw,color = color)

	## PLOT NUMBER OF SEAM CELLS
	if plotNSeam:
	    # create right axes
	    axes2 = axes.twinx()

	    times = []
	    scn = []
	    for t in df.times:
	        if len(df.ix[df.times==t]) > 1:
	            times.append(t)
	            scn.append(len(df.ix[df.times==t]) - 1)
	    pl2, = axes2.plot( times, scn, 'b-', lw = 2 )

	    for tl in axes2.get_yticklabels():
	        tl.set_fontsize(15)    

	    axes2.set_ylim(0,30*5)
	    axes2.plot([0,45],[30,30],'--',color = pl2.get_color(), lw = 2)
	    axes2.set_ylabel('# seam cells', fontsize = 15)
	    axes2.set_yticks(np.arange(0,30+1,10))    
	    axes2.yaxis.label.set_color(pl2.get_color())
	    axes2.tick_params(axis='y', colors=pl2.get_color())

	## PLOT LETH DATA
	if plotLeth:
	    ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )

	    index = np.where( ecd[:,0] == float(worm[1:]) )
	    lethtp = ecd[index, 2:][0][0] - np.min(ecd[index, 1:][0][0]) - 1

	    for tidx in lethtp:
	        tp = df.ix[ (df.rowtype=='body') & (df.tidx==tidx), 'times' ].values[0]
	        axes.plot([tp,tp],[0,21],'--', color = color, lw=2)
	# gc.collect()




if __name__ == '__main__':

	### SETUP THE FIGURE
	fig = plt.figure(figsize = (5.5,4.3))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax.get_xticklabels():
	    tl.set_fontsize(18)
	for tl in ax.get_yticklabels():
	    tl.set_fontsize(18)     
	ax.set_ylim([0,16])
	ax.set_yticks([5,10,15])
	ax.set_xlim([0,50])
	ax.set_xticks([10,20,30,40])
	ax.set_xlabel('Hours after hatching', fontsize = 18)
	ax.set_ylabel('Nuclear fluorescence [a.u.]', fontsize = 18)

	################################################################## DEFINE ALL THE WORMS
	paths = [ '', '' ]
	paths[0] = 'Y:\\Images\\150914_mlt10_250x250x20\\'
	paths[1] = 'Y:\\Images\\150824_JVZ20_wrt2GFP_250x250x20\\' 
	worms = [ [], [] ]
	worms[0] = ['C01','C03','C04','C05','C07','C09','C10','C13','C14','C15','C16','C17','C18','C19','C20']
	worms[1] = ['C01','C02','C03','C04','C05','C06','C07','C09','C10','C11','C12','C13','C14','C18','C19','C22','C23','C24','C25','C26','C27','C29','C30']

	totLength = np.sum( [ len(i) for i in worms ] )
	colors = np.array( [cm.jet(int(i))[:3] for i in np.arange(totLength)*254./(totLength-1)] )
	colors = [ [ i for i in colors[:len(worms[0])] ], [ i for i in colors[len(worms[0]):]] ]

	correctForStrain = False
	correctForSingleWorms = False

	############################################################### PLOT mlt10 DATA
	gene = 'mlt10'
	path = paths[0]
	ws = worms[0]
	colors = colors[0]
	
	# ALL WORMS
	for idx, worm in enumerate( ws ):
	    print( path, worm )
	    plotExpressionOneWorm( ax, gene, path, worm, correctForStrain , correctForSingleWorms,
	    					plotDots=True, plotMean=False, plotSmooth=False, plotLeth = False, plotNSeam = False, 
	    					color = colors[idx], colorDots = colors[idx], alpha = .8)
	
	# # SINGLE WORM DATA
	# worm = 'C10'
	# print( path, worm )
	# plotExpressionOneWorm( ax, 'mlt10', path, worm, correctForStrain, correctForSingleWorms, 
	# 						plotDots=True, plotMean=True, plotSmooth=False, plotLeth = True, plotNSeam = False, 
	# 						color = [0,0,0,1], alpha = 1.)

	# # ############################################################### PLOT wrt2 DATA
	# gene = 'wrt2'
	# path = paths[1]
	# ws = worms[1]
	# colors = colors[1]
	
	# # # ALL WORMS
	# # for idx, worm in enumerate( ws ):
	# #     print( path, worm )
	# #     plotExpressionOneWorm( ax, gene, path, worm, correctForStrain, correctForSingleWorms,
	# # 						plotDots=False, plotMean=False, plotSmooth=True, plotLeth = True, plotNSeam = False, 
	# # 						color = colors[1][idx] )
	
	# # SINGLE WORM DATA
	# worm = 'C05'
	# print( path, worm )
	# plotExpressionOneWorm( ax, 'wrt2', path, worm, correctForStrain, correctForSingleWorms, 
	# 						plotDots=True, plotMean=False, plotSmooth=True, plotLeth = True, plotNSeam = False, 
	# 						color = [0,0,0,1], alpha = 1.)

	# # # # SINGLE CELL MULTIPLE WORMS
	# # ws = ws[::4]
	# # colors = ['k','r','g','b','c','m']
	# # cell = '1.'
	# # shift = 15*np.arange(len(ws))
	# # for idx, worm in enumerate( ws ):
	# #     print( path, worm )
	# #     plotExpressionOneWorm( ax, gene, path, worm, correctForStrain, correctForSingleWorms,
	# # 						plotDots=True, plotMean=False, plotSmooth=False, plotLeth = False, plotNSeam = False, 
	# # 						color = colors[idx], lineage = cell, shift = shift[idx] )

	plt.show()
