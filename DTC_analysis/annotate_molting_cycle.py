import tifffile as tif
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import os.path

indir='W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C10\\'


channel='CoolLED'

# keys to save (h)atch, L(1-4) ecdysis
key=['h', '1', '2', '3', '4']

# keys to move in time
key_time=['q','w','e','r']
DT=[-1,1,-10,10]

class dialog_window:
    def __init__(self, fig):
        self.fig = fig
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.t=0
        self.show_image()
        
    def __call__(self, event):
        # q - increase t, w - decrease t
        # e - increase t by 10, r - decrease t by 10
        if event.key in key_time:
            # change t
            ind=key_time.index(event.key)
            self.t=self.t+DT[ind]
            # check bounds are respected
            if self.t<0:
                self.t=0
            elif self.t>=im.shape[0]:
                self.t=im.shape[0]-1
            # show image
            self.show_image()
        elif event.key in key:
            # save time of (h)atch or L(1-4) ecdysis
            ind=key.index(event.key)
            T[ind]=self.t      
            # show image
            self.show_image()
        elif event.key == ' ':
            # add new event in order (h)atch, L(1-4) ecdysis
            cont=True
            ind=0
            while cont:
                # search for first non-modified event
                if T[ind]==-1:
                    # if found, stop search and assign frame number
                    cont=False
                    T[ind]=self.t
                ind+=1
                # if no unmodified event, stop search
                if ind>=len(T):
                    cont=False
            # save image
            self.show_image()
        elif event.key == 'z' or event.key == 'x':
            # z/x - find closest event before/after current frame
            if event.key=='z':
                # search for earlier
                dt=-1
            else:
                # search for later
                dt=1
            cont=True
            t=self.t
            while cont:
                t=t+dt
                # check if frame already defined
                if t in T:
                    # if so, stop search and move current frame to this frame
                    cont=False
                    self.t=t
                elif t<0 or t>im.shape[0]:
                    # if nothing defined in this direction, stop search
                    cont=False
            self.show_image()
        elif event.key=='escape':
#            # save data
            outfile = indir + "analyzed_data\\molting_cycle_data.p"
            pickle.dump( T, open(outfile, "wb" ) )
            # close figure and exit
            plt.close(self.fig)
            
        self.fig.canvas.draw()
        
    def show_image(self):
        # clear figure
        self.fig.clf()
        W=im.shape[1]
        H=im.shape[2]
        # image fills entire window, no axes
        plt.axes([0,0,1,1])
        # show image
        plt.imshow(im[self.t],cmap='gray', interpolation='nearest')
        # indicate frame
        plt.text(10,30,"%d"% self.t, color='r', fontsize=36)
        # draw a square around each frame that is a hatch or ecdysis
        for i in range(0,len(T)):
            if T[i]==self.t:
                plt.plot([1,W-1,W-1,1,1],[1,1,H-1,H-1,1],'-r',lw=5)
        # print frames of hatch and ecdyses 
        s=""
        for i in range(0,len(T)):
            if T[i] != -1:
                s=s+key[i]+":%d " % T[i]
        if len(s)>0:
            plt.text(10,H-10,s,color='r',fontsize=14)
        # set axis limits
        plt.xlim([0,W])
        plt.ylim([H,0])

# read downsampled trans timelapse stack
filename = indir + "analyzed_data\\"+ channel + "_small.tif"
with tif.TiffFile(filename) as f:
    im=f.asarray()

# load existing data, if there
outfile = indir + "analyzed_data\\molting_cycle_data.p"
if os.path.isfile(outfile):
    T=pickle.load( open(outfile, "rb" ) )
else:
    # else set everything to -1
    T=np.array([-1 for x in range(0,len(key))])


# make figure
fig=plt.figure(num=1, figsize=(6, 6))
# get data
tmp=dialog_window(fig)


