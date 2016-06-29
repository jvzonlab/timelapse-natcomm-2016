import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import tifffile as tif
import re
import os

from POI_fun import POI, find_animal_center
import spline_projection as spp

class animal:
    def __init__(self,x,st,z0,ax,L,lbl):
        self.x=x
        self.st=st
        self.z0=z0
        self.ax=ax
        self.L=L
        self.lbl=lbl

data_set=['W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C01',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C02',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C03',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C05',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C10',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C11',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C12',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C13',
       'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C14']

HW=500
HH=50
f_sc=2

channel='488nm'
track_cell_lbl='d'

for indir in data_set:

    tmp = os.listdir(indir+'\\analyzed_data\\')
    n_list=[]
    for i in range(0,len(tmp)):
        match=re.search('POI_(\d\d\d).p', tmp[i])
        if match:
            n_list.append(int(match.group(1)))
    
    # sort list
    n_list.sort()
    
    # get the information on POIs for all worms in a single list
    animal_data=[]
    for n in n_list:
        # load POIs
        infile = indir + "\\analyzed_data\\POI_%03d.p" % n 
        POI_list=pickle.load( open(infile, "rb" ) )
        
        # find central slice containing POIs 1,2/v and 3
        z0=find_animal_center(POI_list)
    
        infile = indir + "\\analyzed_data\\axis_%03d.p" % n 
        # load axis data, if exists
        ax=pickle.load( open(infile, "rb" ) )
        
        # get spline
        (L,spline_x,spline_y)=spp.get_axis_spline(ax)
        
        # for all cells, get label, (x,y,z) and (s,t) positions
        lbl=[]
        x=[]
        st=[]
        for p in POI_list:
            lbl.append(p.lbl)
            x.append(p.x)
            st.append(spp.get_spline_coordinates_st_for_x(p.x[0:2],spline_x,spline_y,L))
        animal_data.append(animal(np.array(x),np.array(st),z0,ax,L,lbl))
    
    for i in range(0,len(n_list)):
        n=n_list[i]
        
        outfile = indir + "\\analyzed_data\\" + channel + "_straight_%03d.tif" % n
        print "Creating %s" % outfile
    
        # find z slices containing the cells that are tracked
        track_cell_z=[]
        for j in range(0,len(animal_data[i].lbl)):
            if animal_data[i].lbl[j] in track_cell_lbl:
                track_cell_z.append(int(animal_data[i].x[j,2]))
    
        # for now, define (t0,s0), corresponding to the center of the straightened image
        # so that body axis and the half-worm position are at the center
        t0=0
        s0=animal_data[i].L/2
    
        # find t coordinate for gonad, this is needed to know left-right orientation
        for j in range(0,len(animal_data[i].lbl)):
            if int(animal_data[i].x[j,2])==animal_data[i].z0 and animal_data[i].lbl[j] in ['2','v']:
                t_gonad = animal_data[i].st[j,1]
                if animal_data[i].lbl[j] == 'v':
                    # if the vulva is visible, then use its center as the reference point for DTC position
                    s0=animal_data[i].st[j,0]
    
        # for each z slice that will be shown, check whether there are POIs that are also
        # defined in the slize z0. If so use their s coordinate to calculate an offset ds0
        ds0=np.zeros(len(track_cell_z))
        # for each z slice to be shown
        for u in range(0,len(track_cell_z)):
            if track_cell_z[u] != animal_data[i].z0:
                # if not in central slice (where s0 is determined), loop over all POIs
                for v in range(0,len(animal_data[i].lbl)):
                    # find those POIS in current slice that are not tracked
                    if animal_data[i].x[v,2]==track_cell_z[u] and animal_data[i].lbl[v] not in track_cell_lbl:
                        for w in range(0,len(animal_data[i].lbl)):
                            # then, find the point with the same label in slice z0
                            if animal_data[i].lbl[w]==animal_data[i].lbl[v] and animal_data[i].x[w,2]==animal_data[i].z0:
                                # and calculate the difference in s. this is the offset
                                ds0[u]=animal_data[i].st[w,0]-animal_data[i].st[v,0]
    
        # get animals body axis as spline    
        (L,spline_x,spline_y)=spp.get_axis_spline(animal_data[i].ax)
        
        # read in images stack
        im_tracked=np.zeros((2048,2048,len(track_cell_z)))
        filename = indir + "\\z%03d_" % n + channel +'.tif'
        # read in image for each cell tracked
        for j in range(0,len(track_cell_z)):
            im_tracked[:,:,j]=tif.imread(filename, key=track_cell_z[j])    
            
        tmp=np.zeros((2*HH+1,2*HW+1,len(track_cell_z)))
        for j in range(0,len(track_cell_z)):
            # loop over all x in straightened image
            for x in range(0,2*HW+1):
                # for each x, calculate s in (s,t)-projection
                s=s0-ds0[j]+f_sc*(x-HW-1)
                # if s is inside the animal (0<=s<L), then calculate x_s position on spline
                if (s>=0) & (s<=L):
                    # calculate (x(s),y(s)) of point on spline for current value of s
                    x_s=np.array([spline_x(s), spline_y(s)])
                    # get unit vectors e_s and e_t in local (s,t)-coordinate system
                    e_s=np.array([spline_y.derivatives(s)[1],spline_x.derivatives(s)[1]])
                    e_s=e_s/np.sqrt(np.sum(e_s**2))
                    e_t=np.array([e_s[1], -e_s[0]])
                    
                    # loop over all y in straightened image
                    for y in range(0,2*HH+1):
                        # for each y, calculate t in (s,t)-projection
                        t=t0+f_sc*(y-HH-1)
                        # now, calculate coordinates (X,Y) in original image
                        X=int(x_s[0]+t*e_t[1])
                        Y=int(x_s[1]+t*e_t[0])
                        if (X>=0) & (X<im_tracked.shape[1]) & (Y>=0) & (Y<im_tracked.shape[0]):
                            tmp[y,x,j]=im_tracked[Y,X,j]
    
        # get maximum projection of all z slices shown                         
        im_straight=np.max(tmp,axis=2)
        
        # now, depending on the left-right orientation of animal flip image if necessary
        if t_gonad>0:
            im_straight=np.flipud(im_straight)
        
        tif.imsave(outfile,np.array(im_straight, dtype='uint16'))