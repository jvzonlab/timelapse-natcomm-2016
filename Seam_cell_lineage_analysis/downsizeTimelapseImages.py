# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 15:36:52 2015

@author: gritti
"""

import os
import shutil
import glob
import numpy as np
import time
from timelapseFun import *

####################################################################################################
# Give the path and the worm name to be converted
####################################################################################################

path = 'X:\\Nicola\\N2_on_plates\\'

wnumber = np.arange(1,17)
worms = ['C01']
print(worms)

scaleFactor = 8

if type( worms ) == str:
    worms = list(worms)

####################################################################################################
# Downsize function
####################################################################################################

def downsize():
    for worm in worms:
        
        print('Downsizing worm '+worm)

        ################################################################################################
        # Load the list of the images and metadata files
        ################################################################################################
        
        # define the input folder and the output folder
        inpath = path+worm
        outpath = path+worm+'_downsized'
        
        # read in all the images to be converted
        channels = [False,False,False]
        flist = [None,None,None]
        
        if os.path.isfile(inpath+'\\z001_488nm.tif'):
            flist[0] = glob.glob(inpath+'\\z*488nm.tif')
            flist[0].sort()
            channels[0] = True
        if os.path.isfile(inpath+'\\z001_561nm.tif'):
            flist[1] = glob.glob(inpath+'\\z*561nm.tif')
            flist[1].sort()
            channels[1] = True
        if os.path.isfile(inpath+'\\z001_CoolLED.tif'):
            flist[2] = glob.glob(inpath+'\\z*CoolLED.tif')
            flist[2].sort()
            channels[2] = True
        
        metalist = glob.glob(inpath+'\\z*.txt')
        metalist.sort()
        
        ################################################################################################
        # Create the directory and copy metadata files
        ################################################################################################
        
        # create the directory
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        
        # copy the metadataFiles
        for f in metalist:
            # print('copying ' + f)
            if not os.path.isfile(outpath+'\\'+f.split('\\')[-1]):
                shutil.copyfile(f, outpath+'\\'+f.split('\\')[-1])
        
        ################################################################################################
        # Downsize and save the images
        ################################################################################################
        
        # for each channel, if True
        for jdx, chn in enumerate( channels ):
        
            if chn:
                
                # for each image in that channel
                for idx in np.arange(len(flist[jdx])):
                    
                    
                    f = flist[jdx][idx]
                    if not os.path.isfile(outpath+'\\'+f.split('\\')[-1]):
                
                        print('loading ' + flist[jdx][idx])
                        
                        stack = loadstack(f)
                        
                        if len(stack.shape) == 3:
                            
                            smallstack = []
                            for img in stack:
                                Nbig = img.shape[0]
                                Nsmall = img.shape[0]/scaleFactor
                                smallimg = ( img.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1) ).astype(np.uint16)
                                smallstack.append( smallimg )
                            smallstack = np.array(smallstack)

                        if len(stack.shape) == 2:

                            Nbig = stack.shape[0]
                            Nsmall = stack.shape[0]/scaleFactor
                            smallstack = ( stack.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1) ).astype(np.uint16)

                        smallstack = np.array(smallstack)
                    
                        imsave(outpath+'\\'+f.split('\\')[-1],smallstack)
        
    print('Waiting for one hour to restart!')
    time.sleep(60*60)

####################################################################################################
# KEEP DOWNSIZING IMAGES UNTIL SCRIPT IS INTERRUPTED
####################################################################################################

while True:
    downsize()