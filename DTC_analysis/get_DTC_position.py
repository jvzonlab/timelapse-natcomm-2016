import numpy as np
import cPickle as pickle
import re
import os
import datetime as dt

from POI import POI, find_animal_center
import spline_projection as spp

def get_datetime_for_timepoint(indir, n):
    # construct filename
    infile=indir + "z%03d_" % n + '.txt'
    # read text
    text_file = open(infile, "r")
    lines = text_file.readlines()
    # find line with datetime info
    datetime_found=False
    i=0
    while datetime_found==False and i<len(lines):
        # if matches heading of datetime line
        match=re.search('Date/Time: (\S*\s\S*)', lines[i])
        if match:
            # get the datetime and stop looking
            datetime_str=match.group(1)
            datetime_found=True
        i+=1
    # return datetime info
    return (dt.datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S"))

class cell_data:
    def __init__(self,N_T,N_id):
        # N_T - number of timepoints
        # N_id - number of cells tracked in (s,t)-coordinates        
        
        # time since L1 ecdysis
        self.T=np.zeros(N_T,dtype=float)
        # stack number corresponding to each time point
        self.n=np.zeros(N_T,dtype=int)
        # length of worm (anus-posterior pharynx) at each time point
        self.L=np.zeros(N_T,dtype=int)
        # (s,t)
        self.s=np.zeros((N_T,N_id))
        self.t=np.zeros((N_T,N_id))
        
def get_cell_data_for_expt(indir,correct_for_movement):
    
    tmp = os.listdir(indir+'analyzed_data\\')
    n_list=[]
    for i in range(0,len(tmp)):
        match=re.search('POI_(\d\d\d).p', tmp[i])
        if match:
            n_list.append(int(match.group(1)))
    
    # find time point 0, time of L1 ecdysis
    in_file=indir +'analyzed_data\\molting_cycle_data.p'
    tmp=pickle.load( open(in_file, "rb" ) )
    T0=tmp[1]+1
    datetime0=get_datetime_for_timepoint(indir,T0)
    
    DTC_data=cell_data(len(n_list),2)
    
    nn=0
    for n in n_list:
        # get time relative to L1 ecdysis in minutes
        # first, get timestamp
        datetime=get_datetime_for_timepoint(indir,n)
        # calculate difference
        diff=datetime-datetime0
        # convert to minutes
        T=(diff.days * 24 * 60) + (diff.seconds/60)
    
        # load POIs
        outfile = indir + "\\analyzed_data\\POI_%03d.p" % n 
        POI_list=pickle.load( open(outfile, "rb" ) )
        
        # find central slice containing POIs 1,2/v and 3
        Z0=find_animal_center(POI_list)
        
        infile = indir + "\\analyzed_data\\axis_%03d.p" % n 
        print infile
        # load axis data, if exists
        ax=pickle.load( open(infile, "rb" ) )
        
        # get spline
        (L,spline_x,spline_y)=spp.get_axis_spline(ax)
        
        # find t coordinate for gonad
        vulva_midpoint_exist=False
        for i in range(0,len(POI_list)):
            if int(POI_list[i].x[2])==Z0 and POI_list[i].lbl in ['2','v']:
                (s_gonad,t_gonad)=spp.get_spline_coordinates_st_for_x(POI_list[i].x[0:2],spline_x,spline_y,L)
                if POI_list[i].lbl == 'v':
                    # if the vulva is visible, then use its center as the reference point for DTC position
                    vulva_midpoint_exist=True
                    
        # find DTC coordinates (s,t) in spline coordinate system
        tmp=np.zeros((2,4))
        c=0
        for i in range(0,len(POI_list)):
            if POI_list[i].lbl == 'd':
                # get slice number Z
                Z=int(POI_list[i].x[2])
                # get (s,t)-coordinates
                (s,t)=spp.get_spline_coordinates_st_for_x(POI_list[i].x[0:2],spline_x,spline_y,L)
                # by our definition, the position of the gonad should have t<0. If not, then t -> -t for all POIs
                if t_gonad>0:
                    t=-t
                    z=Z-Z0
                else:
                    z=Z0-Z
                # add to array
                tmp[c,:]=np.array([s,t,z,Z])
                c+=1
        # now sort DTC position based on the value of the s coordinate, this will ensure
        # that the first element is the anterior DTC and the second the posterior DTC
        tmp=tmp[np.argsort(tmp[:, 0])]
        
        # find center of animal
        if vulva_midpoint_exist:
            # if the vulva is there, use that as the center (POI 'v')
            s0=s_gonad
        else:
            # else, just use half the measured body length
            s0=L/2
        
        if correct_for_movement==True:                   
            for i in range(0,len(POI_list)):
                # find 'double'anatomical markers. If a POI is indicated in multiple slices
                # that is because the animal moved between these slices
                if POI_list[i].lbl in ['1','2','3','v'] and int(POI_list[i].x[2])!=Z0:
                    lbl_mark=POI_list[i].lbl
                    Z_mark=int(POI_list[i].x[2])
                    (s_mark,t_mark)=spp.get_spline_coordinates_st_for_x(POI_list[i].x[0:2],spline_x,spline_y,L)
                    # find 'original' POI in slice Z0
                    for j in range(0,len(POI_list)):
                        if POI_list[j].lbl==lbl_mark and int(POI_list[j].x[2])==Z0:
                            (s_mark_0,t_mark_0)=spp.get_spline_coordinates_st_for_x(POI_list[j].x[0:2],spline_x,spline_y,L)
                            ds = s_mark - s_mark_0
                    # find DTC in slice Z_mark
                    for j in range(0,2):
                        if int(tmp[j,3])==Z_mark:
                            tmp[j,0]+=ds
    
        DTC_data.s[nn,:]=tmp[:,0]-s0
        DTC_data.t[nn,:]=tmp[:,1]
        DTC_data.T[nn]=T
        DTC_data.L[nn]=L
        DTC_data.n[nn]=n
        nn=nn+1
        
    # sort data in time
    ind_sort=np.argsort(DTC_data.T)
    DTC_data.T=DTC_data.T[ind_sort]
    DTC_data.L=DTC_data.L[ind_sort]
    DTC_data.n=DTC_data.n[ind_sort]
    DTC_data.s=DTC_data.s[ind_sort,:]
    DTC_data.t=DTC_data.t[ind_sort,:]

    return DTC_data

def get_DTC_position(data_dir_list,correct_for_movement):
    DTC_data=[]
    for s in data_dir_list:
        indir=s+'\\'
    
        print "Reading data from %s" % indir
        DTC_data.append(get_cell_data_for_expt(indir,correct_for_movement))
    return DTC_data
