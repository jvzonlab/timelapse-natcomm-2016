# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:00:56 2015

@author: gritti
"""

import glob
from timelapseFun import *
from tifffile import *
import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def isDaughter(cname, cdata):
    return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() == (len(cname)+1) ) ).any()


# path = 'Y:\\Images\\151124_EMS3_250x250x20\\'
# worms = ['C19']

path = 'Y:\\Images\\150325_seamCell_JR667\\'
worms = ['C13']

# path = 'Y:\\Images\\150824_JVZ20_wrt2GFP_250x250x20\\'
# worms = ['C01']

# worms = ['C01','C02','C03','C04','C05','C06','C07','C09','C10','C12','C13','C14','C18','C19','C25','C26','C27','C29','C30']
mothercells = ['1.','2.','3.','4.']#['1.','2.','3.','4.','5.']
cside = 'L'

## for scm data
marginUp = 50
marginDown = 50 
marginLeft = 550
marginRight = 50
marginZ = 1

# ## for wrt-2 data
# marginUp = 25
# marginDown = 25 
# marginLeft = 50
# marginRight = 50
# marginZ = 1

timestep = 1

for w in worms:
    # w = 'C01'
    # outpath = 'Y:\\Presentation\\Figures\\151130_seamcells\\movie_150617_JVZ20_wrt2GFP_%s' % w
    outpath = 'Y:\\Presentation\\Figures\\151130_seamcells\\movie_' + path.split('\\')[-2] + '_%s' %w
    print(path.split('\\'),outpath)

    # print(outpath)
    
    if not os.path.exists(outpath):
        os.mkdir( outpath )
    
    with open(path+'\\worm'+w+'.pickle','rb') as fo:
        df = pickle.load(fo)
        # cdata = df.ix[ (df.rowtype=='cell')&(df.cside==cside), ['cname','times','cside'] ]
        cdata = df.ix[ (df.rowtype=='cell'), ['cname','times','cside'] ]
    print(cdata)
    finalSeamCells = cdata.ix[cdata.times==np.max(cdata.times.values),'cname'].values

    flist = glob.glob(path + w + '_straighten\\straight*_488nm.tif')
    flist.sort()
    flistnew = ['' for i in df.ix[(df.rowtype=='body')&(df.tidx>=0),'tidx']]

    firsttidx = int(flist[0].split('_')[-2][-3:]) - 1
    for f in flist:
        tidx = int(f.split('_')[-2][-3:])-1
        flistnew[tidx-firsttidx] = f

    for i in np.arange(np.min(df.ix[df.rowtype=='cell','tidx'])):
        flistnew.insert(0,'')

    cellrowmask = ( df.rowtype == 'cell' )
    csidemask = df.cside == cside

    bodydata = df.ix[ df.rowtype == 'body' ]
    
    # print(bodydata.ix[bodydata.tidx==0])

    movie = [ [] for i in mothercells ]
    for tidx in bodydata.ix[ bodydata.tidx>=0, 'tidx' ][::timestep]:
        tpmask = df.tidx == tidx
        print(tidx,flistnew[int(tidx)],len(df.ix[tpmask&cellrowmask]))

        for idx, mc in enumerate( mothercells ):
            mothercellmask = df.cname.str.startswith( mc )
            cells = df.ix[ tpmask & mothercellmask & csidemask ].reset_index(True)
 
            if len( cells > 0 ):
                cell = cells.ix[ cells.index[ len( cells.index ) - 1 ] ]
                # print(cell)

                img = loadstack(flistnew[int(tidx)])

                smallimg = img[ cell.cZpos,
                               cell.cYpos - marginUp : cell.cYpos + marginDown,
                               cell.cXpos - marginLeft : cell.cXpos + marginRight ]

                # # max projection
                # Zmin = np.min(cells.cZpos)
                # Zmax = np.max(cells.cZpos)
                # smallimg = np.max( img[ Zmin:Zmax+1,
                #               cell.cYpos-marginUp:cell.cYpos+marginDown,
                #               cell.cXpos-marginLeft:cell.cXpos+marginRight ], axis=0)
    
                # position of the arrows -- only for mlt-10 data!!!
                pos = np.array([cell.cXpos,cell.cYpos])

                ar = smallimg.shape[1]/smallimg.shape[0]
                
                fig = plt.figure(figsize = (4.,4./ar))
                ax = fig.add_subplot(111)
    
                imgplot = ax.imshow(smallimg,cmap=plt.gray(), interpolation='none')
                imgplot.set_clim(0,2**14)
                
                ax.autoscale(False)
                ax.axis('Off')
                fig.subplots_adjust(left=0., right=1., top=1., bottom=0.)
                
                # if cell.cname.startswith('1.'):
                ax.text(2,35,'%dh'%int(cell.times),color='white',fontsize=20)

                # if cell.cname.startswith('4.'):
                #     ax.text(2,45,cell.cside,color='white',fontsize=40)

                # only for mlt-10 data!!!
                for jdx, othercell in cells.iterrows():
                    if isDaughter(othercell.cname, cdata) or (othercell.cname in finalSeamCells):
                        ax.annotate('', fontsize=20, xy=(marginLeft - (pos[0]-othercell.cXpos), 80.),
                                   xycoords='data', xytext=(0, -10.),
                                   textcoords='offset points',
                                   arrowprops=dict(width = 1.,
                                                   headwidth = 8.,
                                                   frac = 1.,
                                                   shrink = 0.05,
                                                   linewidth = 1,
                                                   color = 'yellow')
                                   )

                plt.gray()
                fig.savefig(outpath+'\\movieV%s_%s_tp%.3d.tif'%(mothercells[idx][0],cside,cell.tidx))#, dpi = 300)
                plt.gcf()
                plt.close()
                

    # build single cells grayscale movie
    for idx, sc in enumerate( mothercells ):
        print(sc)
        movie = []
        flist = glob.glob(outpath+'\\movieV%s_%s_tp*.tif'%(sc[0],cside))
        flist.sort()
        for jdx, f in enumerate( flist ):
            img = loadstack(f)[:,:,0]
#            os.remove(f)
            movie.append(img)
            
        imsave( outpath + '\\movieV%s_%s_total.tif' % ( sc[0], cside ), np.array( movie ) )

    for mc in mothercells:
        flist = glob.glob(outpath+'\\movieV%s_%s_tp*.tif'%(mc[0],cside))
        flist.sort()
        for f in flist:
            os.remove(f)
    
    # # build the total movie - only for wrt2 data!!!
    # totmovie = [] #[ np.zeros( (1,200) ) for i in np.arange(0,172)]#np.max(df.tidx)) ]
    # movie = []

    # for idx,sc in enumerate( mothercells ):
    #     flist = glob.glob(outpath + '\\movieV%s_%s_total.tif'%(sc[0],cside))[0]
    #     movie.append(loadstack(flist))
    # [ print(i.shape) for i in movie ]
    # for tidx in np.arange(len(movie[0])):

    #     frame = np.concatenate((movie[0][tidx],movie[1][tidx],movie[2][tidx],movie[3][tidx]),axis=0)
    #     print(movie[0][tidx].shape)
    #     print(frame.shape)
    #     totmovie.append(frame)

    # # print(totmovie)
    # totmovie = np.array(totmovie).astype(np.uint8)
    # imsave(outpath+'\\movieSeamCells.tif',totmovie)

plt.show()