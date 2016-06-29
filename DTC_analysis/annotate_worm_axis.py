import tifffile as tif
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import os.path
import re

from scipy import interpolate
from POI import POI, find_animal_center

indir='W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C03\\'
n=111

def calculate_axis_length(ax):
    L=0
    # for each point on axis
    for i in range(0,ax.shape[0]-1):
        # calculate distance to next point
        dx=(ax[i+1,:]-ax[i,:])**2
        # add it to total length
        L+=np.sqrt(dx.sum())
    return L

def add_point_to_axis(ax,x):
    if ax.shape[0]==2:
        # if only the end points (POI 1 and POI 3) are there, insert in the middle
        ax=np.insert(ax,1,x,axis=0)    
    else: 
        # else, insert in the right order, i.e. that which minimizes axis length
        L_min=6e66
        ind_min=-1
        for i in range(1,ax.shape[0]):
            # for each possible point of insertion
            tmp=np.insert(ax,i,x,axis=0)
            # calculate the length
            L=calculate_axis_length(tmp)
            if L<L_min:
                # and save the index of the insert point that results in the minimal total length
                L_min=L
                ind_min=i
        # then insert point at this position
        ax=np.insert(ax,ind_min,x,axis=0)    
    return (ax)

def find_axis_point_closest_to_x(ax,x):
    dx2_min=6e66
    ind_min=-1
    # for each point on axis
    for i in range(1,ax.shape[0]-1):
        # calculate distance^2 to point x
        tmp=(x-ax[i,:])**2
        dx2=tmp.sum()
        if dx2<dx2_min:
            # keep index of closest point so far
            dx2_min=dx2
            ind_min=i
    # and return it
    return (ind_min)

def get_axis_spline(ax):
    # get arc length of axis
    s=np.zeros(ax.shape[0])
    for i in range(1,ax.shape[0]):
        dx=ax[i,:]-ax[i-1,:]
        s[i]=s[i-1]+np.sqrt(np.sum(dx**2))
    # then get spline as function of arc length
    f_spline=1e1
    spline_x = interpolate.UnivariateSpline(s, ax[:, 0],s=f_spline*ax.shape[0])
    spline_y = interpolate.UnivariateSpline(s, ax[:, 1],s=f_spline*ax.shape[0])
    # and return it
    return (s[-1],spline_x,spline_y)
        
class axis_dialog_window:
    def __init__(self, fig, ax):
        self.fig = fig
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.cid_mouse = self.fig.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.ax=ax

        self.show_image()
        
    def on_key_press(self, event):
        if event.key=='escape':
            # save data
            outfile = indir + "\\analyzed_data\\axis_%03d.p" % m 
            pickle.dump( self.ax, open(outfile, "wb" ) )
            # close figure and exit
            plt.close(self.fig)

    def on_button_press(self, event):
        x=np.array([event.xdata*f_sc, event.ydata*f_sc])
        if event.button==1:
            # add point to axis
            self.ax=add_point_to_axis(self.ax,x)
        elif event.button==3:
            # find closest point to mouse cursor
            ind=find_axis_point_closest_to_x(self.ax,x)
            if ind!=-1:
                # if there is a deletable point (i.e. not POI 1 or 3), remove
                self.ax=np.delete(self.ax,ind,axis=0)
        # show result
        self.show_image()

    def show_image(self):
        self.fig.clf()
        # make image fill whole figure
        plt.axes([0,0,1,1])
        # plot trans image
        plt.imshow(data,cmap='gray', interpolation='nearest')
        # plot axis points
        plt.plot(self.ax[:,0]/f_sc,self.ax[:,1]/f_sc, 'or')
        # if sufficient points to fit spline, show it
        if self.ax.shape[0]>3:
            # get the spline fit
            (L,sx,sy)=get_axis_spline(self.ax)
            # then plot it
            s = np.arange(0,L,10)
            plt.plot(sx(s)/f_sc,sy(s)/f_sc, '-r')
        # set limits
        plt.xlim([0,data.shape[0]])
        plt.ylim([0,data.shape[1]])

# get axis in trans image
channel='CoolLED'
# scale image down by factor f_sc to speed up annotation
f_sc=2

n_list=[]   
if 'n' in globals():
    # if n is defined, only analyze animal <n>
    n_list=[n]
else:
    # else, get first animal for which POI but not axis data exists
    tmp = os.listdir(indir+'analyzed_data\\')
    cont=True
    i=0
    while cont and i<len(tmp):
        match=re.search('POI_(\d\d\d).p', tmp[i])
        if match:
            m=int(match.group(1))
            if os.path.isfile(indir + "\\analyzed_data\\axis_%03d.p" % m )==False:
                n_list=[m]
                print n_list
                cont=False
        i+=1

for m in n_list:     
    # load POIs
    outfile = indir + "\\analyzed_data\\POI_%03d.p" % m 
    print outfile
    POI_list=pickle.load( open(outfile, "rb" ) )
    
    # find central slice containing POIs 1,2/v and 3
    Z0=find_animal_center(POI_list)

    # if not there, quit    
    if Z0==-1:
        raise ValueError("Z0 not properly defined, exiting program.")
        
    outfile = indir + "\\analyzed_data\\axis_%03d.p" % m 
    if os.path.isfile(outfile):
        # load axis data, if exists
        ax=pickle.load( open(outfile, "rb" ) )
    else:
        # else, intialize worm axis: set POIs '1' and '3' at Z0 as start and end of axis
        ax=np.zeros((2,2),dtype=float)
        for i in range(0,len(POI_list)):
            if int(POI_list[i].x[2])==Z0:
                if POI_list[i].lbl=='1':
                    ax[0,:]=POI_list[i].x[0:2]
                elif POI_list[i].lbl=='3':
                    ax[1,:]=POI_list[i].x[0:2]
    
    # read in trans image
    filename = indir + "z%03d_" % m + channel +'.tif'
    data=tif.imread(filename,key=Z0)
    data=data[::f_sc,::f_sc]
    
    # get axis
    fig=plt.figure(num=1, figsize=(11,11))
    tmp=axis_dialog_window(fig,ax)