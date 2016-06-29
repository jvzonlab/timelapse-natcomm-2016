import tifffile as tif
import numpy as np
import os
import re

# use trans channel
channel='CoolLED'
# scale images down by factor 6
f_sc=6

data_dir='W:\\Jeroen\\timelapse_data\\151016_unc6_qIs56\\'

# find all animals in data directory
tmp=os.listdir(data_dir)
for s in tmp:
    # for each animal
    if os.path.isdir(data_dir + s):
        # find all trans images
        animal_dir=data_dir + s + '\\'
    
        # find all trans images
        tmp2 = os.listdir(animal_dir)
        n_list=[]
        for i in range(0,len(tmp2)):
            # that match format z[???]_<channel>.tif
            match=re.search('z(\d\d\d)_' + channel + '.tif', tmp2[i])
            if match:
                n_list.append(int(match.group(1)))
        # sort images in ascending order
        n_list.sort()

        # load all images
        c=0
        for n in n_list:
            filename = animal_dir + "z%03d_" % n + channel +'.tif'
            print "Reading data from: %s" % filename
            
            if c==0:
                # read in trans stack, find number of slices
                with tif.TiffFile(filename) as f:
                    tmp=f.asarray()
                    # sample down <f_sc>-fold
                    tmp=tmp[:,::f_sc,::f_sc]
                # intialize array for stack                    
                data=np.zeros((len(n_list),tmp.shape[1],tmp.shape[2]),dtype='uint16')
                # get center of stack, z=Z
                Z=int(tmp.shape[0]/2)
                # save image
                data[c,:,:]=tmp[Z,:,:]
            else:
                # only read in slice at z=Z
                tmp=tif.imread(filename,key=Z)
                # sample down <f_sc>-fold
                tmp=tmp[::f_sc,::f_sc]
                # save image
                data[c,:,:]=tmp
            c=c+1
    # check if directory exists,
    outdir=animal_dir+"analyzed_data\\"
    if not os.path.exists(outdir):
        # if not, create
        os.makedirs(outdir)            
    # then save stack
    with tif.TiffWriter(outdir + channel+'_small.tif') as f:
        for i in range(data.shape[0]):
            f.save(data[i])
            