import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import os.path

from POI import POI, find_animal_center

N=100

indir='W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C02\\'
N_range=[48,123]
N_skip=[]

# name of trans and GFP channel
channel=['CoolLED', '488nm']
# scale images by this factor for speed 
f_sc=2

def load_stack(n,indir,channel):
    tmp=[]
    data=[]
    for i in range(0,len(channel)):
        # get filename for this time point
        filename = indir + "z%03d_" % n + channel[i] +'.tif'
        print "Reading data from: %s" % filename
    
        # read in trans stack
        with tiff.TiffFile(filename) as f:
            tmp=f.asarray()
            # scale down image
            tmp=tmp[:,::f_sc,::f_sc]
        if i==0:
            # if first channel, then initialize array to store all stacks
            data=np.zeros((len(channel),tmp.shape[0],tmp.shape[1],tmp.shape[2]),dtype=int)
        # save data to array
        data[i,:,:,:]=tmp
        tmp=[]
    return(data)

# keys:
# q - one slice down
# w - one slice up
# e - five slices down
# r - five slices up
# z - view trans channel
# x - view GFP channel
# 1,2,3,d,v - add label at current mouse position and current slice

class POI_dialog_window:
    def __init__(self, fig, POI_list):
        # keys to move in z
        self.key_z=['q','w','e','r']
        self.dz=[-1,1,-5,5]
        # keys to change channel
        self.key_ch=['z','x']
        # POI labels: 1,2/v,3,d
        self.key_lbl=['1','2','3','d','v']
        
        self.fig = fig
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)
        
#        # find min and max intensity for all channels
#        self.int_minmax=np.zeros((data.shape[0],2))
#        for i in range(0,data.shape[0]):
#            self.int_minmax[i,0]=np.min(data[i,:,:,:])
#            self.int_minmax[i,1]=np.max(data[i,:,:,:])
        self.maxint_scale=np.zeros(data.shape[0])
        for i in range(0,data.shape[0]):
            self.maxint_scale[i]=0.7
        
        # initialize channel to trans
        self.ch=0
        # start at center of stack
        self.z=int(data.shape[1]/2)
        # show image
        self.show_image()
        
    def __call__(self, event):
        if event.key in self.key_z:
            # change z
            ind=self.key_z.index(event.key)
            self.z=self.z+self.dz[ind]
            # check bounds are respected
            if self.z<0:
                self.z=0
            elif self.z>=data.shape[1]:
                self.z=data.shape[1]-1
            # show image
            self.show_image()
        elif event.key in self.key_ch:
            # set channel
            self.ch=self.key_ch.index(event.key)
            # show image
            self.show_image()
        elif event.key in self.key_lbl:
            # get label
            lbl=event.key.encode('ASCII')
            # save POI, use coordinates in units of pixels in original image
            POI_list.append(POI([event.xdata*f_sc, event.ydata*f_sc, self.z],lbl))
            # show image
            self.show_image()
        elif event.key=='delete' or event.key=='backspace':
            # remove closest POI
            x=event.xdata
            y=event.ydata
            min_dist_sq=6e66
            min_id=-1
            # for each POI
            for i in range(0,len(POI_list)):
                dist_sq=(x-POI_list[i].x[0]/f_sc)**2+(y-POI_list[i].x[1]/f_sc)**2
                # find the one with the smallest distance to the mouse cursor
                if dist_sq<min_dist_sq:
                    min_dist_sq=dist_sq
                    min_id=i
            # then remove this from the list
            del POI_list[min_id]
            # finally, show image
            self.show_image()
        elif event.key=='[' or event.key==']':
            fac=1.5
            if event.key==']':
                # decrease contrast, i.e. white corresponds to higher pixel value
                self.maxint_scale[self.ch]*=fac
            else:
                # increase contrast, i.e. white corresponds to lower pixel value
                self.maxint_scale[self.ch]/=fac                
            # finally, show image
            self.show_image()
        elif event.key=='escape':
            # find central slice containing POIs 1,2/v and 3
            Z0=find_animal_center(POI_list)
            # save and exit only of 1,2/v and 3 are all defined in a single slice
            if Z0!=-1:
                # save data
                outfile = indir + "\\analyzed_data\\POI_%03d.p" % n 
                pickle.dump( POI_list, open(outfile, "wb" ) )
                # close figure and exit
                plt.close(self.fig)
            else:
                print "ALARMMMM!!!!!"

 
    def show_image(self):
        # clear figure
        self.fig.clf()
        # make axes fill entire figure
        plt.axes([0,0,1,1])
        # calculate color range for image
        int_minmax=(np.min(data[self.ch,self.z,:,:]),self.maxint_scale[self.ch]*np.max(data[self.ch,self.z,:,:]))
        # show image
        plt.imshow(data[self.ch,self.z,:,:],cmap='gray', clim=int_minmax,interpolation='nearest')
        # set limits of axes
        plt.xlim([0,data.shape[2]])
        plt.ylim([0,data.shape[3]])
        # plot all POIs
        for i in range(0,len(POI_list)):
            if POI_list[i].x[2]==self.z:
                # in yellow, if they are in current slice
                col='y'
            else:
                # otherwise, in red
                col='r'
            plt.plot(POI_list[i].x[0]/f_sc,POI_list[i].x[1]/f_sc,'.'+col)
            plt.text(POI_list[i].x[0]/f_sc+10.,POI_list[i].x[1]/f_sc+10.,POI_list[i].lbl,color=col,fontsize=18)

n=-1   
if 'N' in globals():
    # if n is defined, only analyze animal <n>
    n=N
else:
    # else, get first animal for which POI does not exist
    i=N_range[0]
    cont=True
    while cont and i<=N_range[1]:
        tmp=indir + "\\analyzed_data\\POI_%03d.p" % i
        if os.path.isfile(tmp) == False and i not in N_skip:
            n=i
            cont=False
        i+=1

if n>0:        
    # load stack
    data=load_stack(n,indir,channel)
               
    # load POI data, if exists
    outfile = indir + "\\analyzed_data\\POI_%03d.p" % n 
    if os.path.isfile(outfile):
        POI_list=pickle.load( open(outfile, "rb" ) )
    else:
        POI_list=[]
               
    # make dialog window           
    fig=plt.figure(num=1, figsize=(11, 11))
    tmp=POI_dialog_window(fig,POI_list)
