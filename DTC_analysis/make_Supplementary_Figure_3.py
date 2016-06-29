import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import tifffile as tif
import re
import datetime as dt

from get_DTC_position import get_DTC_position

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def get_datetime_for_timepoint(indir,n):
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
    
def get_molting_time(data_set):
    n_molt=np.zeros((len(data_set),5),dtype=int)
    for i in range(0,len(data_set)):
        # for each animal
        in_file=data_set[i] + "\\analyzed_data\\molting_cycle_data.p"
        tmp=pickle.load( open(in_file, "rb" ) )
        n_molt[i,:]=1+tmp
    return (n_molt)

# get sliding average of DTC position x for each timepoint t with window size 2*W
def get_sliding_avg_of_DTC_position(t,x,W):
    x_savg=np.zeros(x.shape)
    # for each time point <t>
    for i in range(0,len(t)):
        # find all data points with t-W <= time <= t+W
        ind=np.where((t>=t[i]-W) & (t<=t[i]+W))
        # get the mean position over all time points
        tmp=np.mean(x[ind])
        # and save
        x_savg[i]=tmp
    return x_savg

# get DTC velocity versus time for all DTCs in <cell_ind>
def get_velocity_vs_time(cell_ind,data,DT_L1L3,T_range,DT_acq,W_savg,axis="DV or AP"):
    # find out how many DTCs in this class of DTCs
    N_cell=len(cell_ind)
    # define timepoints at which to evaluate DTC velocity
    t=np.arange(T_range[0],T_range[1],DT_acq,dtype=float)
    # save array for velocity. -1 corresponds to no value for this cell at this time
    v=-66.*np.ones((len(t),N_cell))
    
    # loop over all DTCs
    for n in range(0,N_cell):
        # get data set
        i=cell_ind[n][0]
        # get animal
        j=cell_ind[n][1]
        # get cell (DTC.a or DTC.p)
        k=cell_ind[n][2]
        if axis=='AP':
            # get sliding average of DTC position s
            x_savg=get_sliding_avg_of_DTC_position(data[i][j].T,data[i][j].s[:,k],W_savg)
        elif axis=='DV':            
            # get sliding average of DTC position t
            x_savg=get_sliding_avg_of_DTC_position(data[i][j].T,data[i][j].t[:,k],W_savg)
        # get time of L3 ecdysis relative to L1 ecdysis
        T0=DT_L1L3[i][j]
        # calculate time relative to L2 ecdysis
        T=data[i][j].T[0:-1]-T0
        # calculate delta_s
        Dx=np.diff(x_savg)
        # calculate delta_time
        DT=np.diff(data[i][j].T)
        # calculate velocity in um per hour
        v_cell=um_per_pixel*Dx/DT
        # for all data points of this DTC, find index in array <t>
        inds=np.array(((T-T_range[0])/DT_acq),dtype=int)
        # save the velocity at the corresponding index in <v>
        if axis=='DV':
            v[inds,n]=v_cell
        if axis=='AP':
            # define velocity so that moving outward corresponds to v>0
            if k==0:
                v[inds,n]=-v_cell
            elif k==1:
                v[inds,n]=v_cell
    return (t,v)

# given the velocity v at time t for all DTCs in a class, 
# calculate the average over the entire class
def get_average_velocity_vs_time(t,v):
    T_avg=[]
    V_avg=[]
    V_stderr=[]
    # for each timepoint
    for i in range(0,v.shape[0]):
        # find all DTC velocities that are defined, i.e. v>-6e66
        ind=np.where(v[i,:]>-66)[0]
        # if any data points found, then
        if len(ind)>0:
            # save the time,
            T_avg.append(t[i])
            # the average velocity
            V_avg.append(np.mean(v[i,ind]))
            # and the standard error
            N=len(ind)
            V_stderr.append(np.std(v[i,ind])/np.sqrt(N-1))
    return ( np.array(T_avg), np.array(V_avg), np.array(V_stderr) )

WT_data_set=['W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C01',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C02',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C03',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C04',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C05',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C06',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C08',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C09',
             'W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\C10']
             
unc6_data_set=['W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C01',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C02',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C03',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C05',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C10',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C11',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C12',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C13',
               'W:\\Jeroen\\timelapse_data\\16_10_2015_unc6_qIs56\\C14']

## load _data
#data=[]
#data.append(get_DTC_position(WT_data_set))
#data.append(get_DTC_position(unc6_data_set))
#pickle.dump(data, open('DTC_data.p', 'wb'))

data = pickle.load(open('figure_data\\DTC_data_movement_corrected.p', 'rb'))

fig=plt.figure(1,figsize=(10,8))
plt.clf()

# get time between L1 molt and L3 molt, so that I can plot DTC dynamics with respect to L3 molt
DT_L1L3=[]
# for each data set
for i in range(0,2):
    # get list of data directories for each animal
    if i==0:
        ds=WT_data_set
    else:
        ds=unc6_data_set
    # get the time of all molts for each data set        
    MT=get_molting_time(ds)
    # calculate the time difference between the L1 and L3 molt in minutes
    DT=[]
    # for each animal
    for j in range(0,len(ds)):
        # get datetime of L1 and L3 molt 
        indir=ds[j]+"\\"
        T_L1=get_datetime_for_timepoint(indir,MT[j,1])
        T_L3=get_datetime_for_timepoint(indir,MT[j,3])
        # calculate difference in minutes and add to list
        diff=T_L3-T_L1
        DT.append((diff.days * 24 * 60) + (diff.seconds/60))
    # add list for each data set    
    DT_L1L3.append(DT)

um_per_pixel=0.1625

# time interval between image acquisitions
DT_acq=20.
# window size for sliding average
W_savg=40
# time range over which to take velocity average
T_range=np.array([-1200,800])

# get indeces of all WT DTCS 
i=0
avg_ind_WT=[]
# for each animal
for j in range(0,len(data[i])):
    # save indeces to both DTCs
    avg_ind_WT.append([i,j,0])
    avg_ind_WT.append([i,j,1])

# get indeces of unc-6 DTCS, split into those that turn around and those that don't
# for those that don't turn, separate in DTC.a and DTC.p 
i=1
avg_ind_unc6_turn=[]
avg_ind_unc6_no_turn_A=[]
avg_ind_unc6_no_turn_P=[]
# for each animal
for j in range(0,len(data[i])):
    # for each DTC
    for k in range(0,2):
        if np.abs(data[i][j].s[-1,k])<(150./um_per_pixel):
            avg_ind_unc6_turn.append([i,j,k]) 
        else:
            print j
            # else, save in "no turn" list...
            if k==0:
                # ...for DTC.a
                avg_ind_unc6_no_turn_A.append([i,j,k])
            else:
                # else, for DTC.p
                avg_ind_unc6_no_turn_P.append([i,j,k])
# collate the index list of all three classes in one list
all_ind_list=[avg_ind_WT,avg_ind_unc6_turn,avg_ind_unc6_no_turn_A,avg_ind_unc6_no_turn_P]

# plot raw and average s vs time date for all three classes
col='rgcm'
plt.subplot(3,3,2)
plt.plot([-15,15],[0,0],'--k')
plt.plot([0,0],[-320,300],'--k')
for m in range(0,len(all_ind_list)):
    ind=all_ind_list[m]
    for n in range(0,len(ind)):
        i=ind[n][0]
        j=ind[n][1]
        k=ind[n][2]
        
        T0=DT_L1L3[i][j]
        plt.plot((data[i][j].T-T0)/60.,um_per_pixel*get_sliding_avg_of_DTC_position(data[i][j].T,data[i][j].s[:,k],W_savg),'-',color=col[m])
plt.xlabel('Hours after L3 ecdysis')
plt.ylabel('A-P position (um)')
plt.xlim([-17,12])
plt.ylim([-320,300])

avg_data=[]
ind_plt=[4,6,7,8]
for m in range(0,len(all_ind_list)):
    (t, v) = get_velocity_vs_time(all_ind_list[m],data,DT_L1L3,T_range,DT_acq,W_savg,axis="AP")
    (t_avg, v_avg, v_stderr) = get_average_velocity_vs_time(t,v)
    avg_data.append([t_avg,v_avg,v_stderr])    
    
    plt.subplot(3,3,ind_plt[m])
    plt.plot(t/60.,v, '.',color=col[m])
    plt.plot(t_avg/60.,v_avg,'-k')
    plt.errorbar(t_avg[::3]/60., v_avg[::3], yerr=v_stderr[::3], fmt='.',color='k')
    plt.plot([0,0],[-1,1],'--k')
    plt.plot(T_range/60.,[0,0],'--k')
    plt.xlim([-17,12])
    plt.ylim([-0.6,0.6])
    plt.xlabel('Hours after L3 ecdysis')
    plt.ylabel('v_A-P (um/hr)')
    
plt.subplot(3,3,9)
plt.plot(T_range/60.,[0,0],'--k')
plt.plot([0,0],[-1,1],'--k')
for m in range(0,len(all_ind_list)):
    t_avg=avg_data[m][0]
    v_avg=avg_data[m][1]
    v_stderr=avg_data[m][2]
    plt.plot(t_avg/60.,v_avg,'-',color=col[m])
    plt.plot(t_avg[::3]/60., v_avg[::3], '.',color=col[m])
#    plt.errorbar(t_avg[::3]/60., v_avg[::3], yerr=v_stderr[::3], fmt='.',color=col[m])
plt.xlim([-17,12])
plt.ylim([-0.3,0.5])
plt.xlabel('Hours after L3 ecdysis')
plt.ylabel('v_A-P (um/hr)')

# plot raw and average t vs time for Wt animals 
plt.subplot(3,3,3)
plt.plot([-15,15],[0,0],'--k')
plt.plot([0,0],[-10,15],'--k')
ind=all_ind_list[0]
for n in range(0,len(ind)):
    i=ind[n][0]
    j=ind[n][1]
    k=ind[n][2]
    
    T0=DT_L1L3[i][j]
    plt.plot((data[i][j].T-T0)/60.,um_per_pixel*get_sliding_avg_of_DTC_position(data[i][j].T,data[i][j].t[:,k],W_savg),'-',color=[0.7,0.7,0.7])
plt.xlabel('Hours after L3 ecdysis')
plt.ylabel('D-V position (um)')
plt.xlim([-17,12])

# calculate dorsal_ventral velocities v_DV
(t, v_DV) = get_velocity_vs_time(all_ind_list[0],data,DT_L1L3,T_range,DT_acq,W_savg,axis="DV")
# get average <v_DV> over all DTCs
(t_avg, v_DV_avg, v_DV_stderr) = get_average_velocity_vs_time(t,v_DV)
# and plot
plt.subplot(3,3,5)
plt.plot(t/60.,v_DV, '.',color=[0.5,0.5,0.5])
plt.plot(t_avg/60.,v_DV_avg,'-k')
plt.errorbar(t_avg[::3]/60., v_DV_avg[::3], yerr=v_DV_stderr[::3], fmt='.',color='k')
plt.plot([0,0],[-1,1],'--k')
plt.plot(T_range/60.,[0,0],'--k')
plt.xlim([-17,12])
plt.ylim([-0.1,0.2])
plt.xlabel('Hours after L3 ecdysis')
plt.ylabel('v_D-V (um/hr)')

# panel A, difference in crossing time
t_c=0
i=0
T_cross=np.zeros((len(data[i]),2))
for j in range(0,len(data[i])):
    for n in range(0,2):    
        tmp=np.where(data[i][j].t[:,n]>t_c)[0]
        c=0
        cont=True
        while cont:
            if c+1<len(tmp):
                if tmp[c+1]==tmp[c]+1:
                    cont=False
                    T_cross[j,n]=data[i][j].T[tmp[c]]
                c+=1
            else:
                cont=False
ax=plt.subplot(331)
plt.plot([0,0],[0,6],'--k')
dt=np.arange(-110,100,20)
h=np.histogram(np.diff(T_cross),dt)
t=h[1][0:-1]
plt.bar(t/60.,h[0],width=0.85*20./60.,color='gray')

plt.ylim([0,6])
plt.xlim([-2,1.5])
plt.setp(ax, xticks=[-2,-1,0,1])
plt.setp(ax, yticks=[])
plt.xlabel('Difference in D-V migration time (hour)')
plt.ylabel('fraction of animals')