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
from matplotlib import cm
mpl.rcParams['pdf.fonttype'] = 42



def plotExpressionOneWorm( axes, path='Y:\\Images\\150617_JVZ20_wrt2GFP\\', worm='C01', plotSingleDots = True, plotLeth = True, plotSmooth = True, color = [0.,0.,0.,1.] ):

	# load the data    
	if os.path.exists( os.getcwd()+'\\worm'+worm+'.pickle' ):
		loadPath = os.getcwd()+'\\worm'+worm+'.pickle'
	else:
		loadPath = path + 'worm' + worm + '.pickle'
	df = pickle.load(open(loadPath,'rb'))
	df.exp /= 1000
	    
	## PLOT THE MEAN
	if plotSingleDots:
		pl1, = axes.plot(df.times[pd.notnull(df.exp)],
		    df.exp[ pd.notnull(df.exp)], 'o', lw=2, mfc = color[:3], alpha = .3 )

	## PLOT THE FILTERED DATA

	if plotSmooth:
		# group data by time and remove missing datapoint
		groups = df.groupby(df.times)
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
				filtdata.filtexp,'-',lw=3,color = color)

	## PLOT LETH DATA
	if plotLeth:
		ecd = np.loadtxt( open( path + 'skin.txt', 'rb' ) )

		index = np.where( ecd[:,0] == float(worm[1:]) )
		lethtp = ecd[index, 2:][0][0] - np.min(ecd[index, 1:][0][0]) - 1

		for tidx in lethtp:
			tp = df.ix[ (df.rowtype=='body') & (df.tidx==tidx), 'times' ].values[0]
			axes.plot([tp,tp],[0,15],'--', color = color, lw = 2)



if __name__ == '__main__':

	### SETUP THE FIGURE
	fig = plt.figure(figsize = (5.8,3.8))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.12, right=.97, top=.97, bottom=0.16)
	for tl in ax.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax.get_yticklabels():
		tl.set_fontsize(18) 
	ax.set_ylim([0,15])
	ax.set_xlim([0,45])
	ax.set_xlabel('Hours after hatching', fontsize = 18)
	ax.set_ylabel('Worm fluorescence [a.u.]', fontsize = 18)

	### DEFINE ALL THE WORMS
	path = 'Y:\\Images\\150914_mlt10_250x250x20\\'
	worms = ['C01','C03','C04','C05','C07','C20']

	colors = np.array( [cm.jet(int(i))[:3] for i in np.arange(len(worms))*254./(len(worms)-1)] )

	# ### PLOT ALL WORMS DATA
	# for idx, worm in enumerate( worms ):
	# 	print( path, worm )
	# 	plotExpressionOneWorm( ax, path, worm, plotSingleDots = False, plotLeth = True, plotSmooth = True, color = colors[idx]  )

	# ### PLOT SINGLE WORM DATA
	worm = 'C05'
	print( path, worm )
	plotExpressionOneWorm( ax, path, worm, plotSingleDots = True, plotLeth = True, plotSmooth = True, color = [0,0,0,1]  )	

	plt.show()