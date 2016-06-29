import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import re
import os
import tifffile as tif

def get_datetime_for_timepoint(n):
    # construct filename
    infile=indir + "z%03d_" % n + '.txt'
    # read text
    text_file = open(infile, "r")
    lines = text_file.readlines()
    # find line with datetime infor
    datetime_found=False
    i=0
    while datetime_found==False and i<len(lines):
        # if matches heading of datetime line
        match=re.search('Date/Time: (\S*\s\S*)', lines[i])
        if match:
            # get the datetike and stop looking
            datetime_str=match.group(1)
            datetime_found=True
        i+=1
    # return datetime info
    return (dt.datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S"))

## nice movies are:    
#'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C01',
#'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C02',
#'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C05',
#'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C11',
#'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C12',

# These are the ones I use:
# wild-type
#indir='W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C10\\'
# both ventral
#indir='W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C11\\'
# one ventral
indir='W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C02\\'

channel='488nm'

# for converting pixels to micrometers
um_per_pixel=0.1625

# get a list of all images of straightened animals
tmp = os.listdir(indir+'analyzed_data\\')
n_list=[]
for i in range(0,len(tmp)):
    match=re.search(channel+'_straight_(\d\d\d)', tmp[i])
    if match:
        n_list.append(int(match.group(1)))
# sort list, so that images occur in ascending order of time point        
n_list.sort()

# get image and time data for each time point
T=np.zeros(len(n_list))
for i in range(0,len(n_list)):
    # read in straightened animal image for current timepoint
    filename=indir + 'analyzed_data\\' + channel + '_straight_%03d.tif' % n_list[i]
    tmp=tif.imread(filename)
    
    if i==0:
        # for first image, get dimensions and initialize array for all time points
        im_data=np.zeros((tmp.shape[0],tmp.shape[1],len(n_list)))    
        # get time for this time point (corresponding to L1 ecdysis)
        datetime0=get_datetime_for_timepoint(n_list[i])
    # save image to array of images
    im_data[:,:,i]=tmp

    datetime=get_datetime_for_timepoint(n_list[i])
    # calculate time difference
    diff=datetime-datetime0
    # and convert to hours
    T[i]=(diff.days * 24.) + (diff.seconds/3600.)

# the GFP intensity increases significantly during time course. In the movie, we want
# to equalize the contrast so that it remains approximately the same for each time point
max_int=np.zeros(len(n_list))
min_int=np.zeros(len(n_list))
# so, for each timepoint
for i in range(0,len(n_list)):
        # find the largest GFP intensity
        max_int[i]=np.max(im_data[:,:,i])
        r=im_data[:,:,0]>0
        # and the smallest non-zero intensity
        min_int[i]=np.min(im_data[r,i])
# then, do a polynomial fit to it to smooth out timepoint-to-timepoint variability        
n=np.array(range(0,len(n_list)))
p=np.polyfit(n, max_int, 4)

my_dpi=100.
for i in range(0,len(n_list)):
    # a pixel value of 0.0 corresponds to the mean lowest GFP intensity
    int_low=np.mean(min_int)
    # pixel value of 1.0 corresponds to 50% of the interpolated max value
    int_high=0.5*np.polyval(p,i)
    
    # clip so that pixel value = 0 for intensities<int_low  
    tmp=im_data[:,:,i]-int_low
    r=tmp<0
    tmp[r]=0
    # then normalize so that int_high corresponds to 1.0
    tmp=tmp/(int_high-int_low)
    
    # plot image
    fig=plt.figure(1,figsize=(tmp.shape[1]/my_dpi, tmp.shape[0]/my_dpi), dpi=my_dpi)
    plt.clf()
    plt.axes([0,0,1,1])
    imgplot=plt.imshow(tmp, cmap='gray', aspect='auto')
    imgplot.set_clim(0.0,1.0)
    # plot annotation, starting with body axis
    plt.plot([0,tmp.shape[1]],((tmp.shape[0]-1)/2)*np.ones(2),':w')
    # plot time
    plt.text(5,30,"%2dh" % np.round(T[i]), color='w', fontsize=24)
    # adjust axes
    plt.xlim([0,tmp.shape[1]])
    plt.ylim([tmp.shape[0],0])
    plt.axis('off')
    
    # this is the scale factor in make_straightened_animal_images
    f_sc=2.
    # plot 20um scale bar
    plt.plot([20,20+20./f_sc/um_per_pixel],[tmp.shape[0]-20,tmp.shape[0]-20],'w',lw=5)

    # save a temporary copy of image as bitmap
    plt.savefig(indir+'analyzed_data\\tmp_%03d.tif' % i, dpi=my_dpi)

# now, load temporary annotated images into numpy array    
for i in range(0,len(n_list)):
    tmp=tif.imread(indir+'analyzed_data\\tmp_%03d.tif' % i)
    if i==0:
        movie_data=np.zeros((len(n_list),tmp.shape[0],tmp.shape[1]),dtype='uint16')
    movie_data[i,:,:]=tmp[:,:,0]

# and save as multipage tif file
tif.imsave(indir+'analyzed_data\\DTC_movie.tif', movie_data)

# then delete the original separate images
for i in range(0,len(n_list)):
    os.remove(indir+'analyzed_data\\tmp_%03d.tif' % i)
      