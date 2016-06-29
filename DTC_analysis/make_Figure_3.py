import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import tifffile as tif
import re
import datetime as dt

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

from POI import POI
from get_DTC_position import get_DTC_position

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
    
def get_molting_frame(data_set):
    n_molt=np.zeros((len(data_set),5),dtype=int)
    for i in range(0,len(data_set)):
        # for each animal
        in_file=data_set[i] + "\\analyzed_data\\molting_cycle_data.p"
        tmp=pickle.load( open(in_file, "rb" ) )
        n_molt[i,:]=1+tmp
    return (n_molt)

def get_ecdysis_time_after_L1_ecdysis(data_set):
    T_ecd=np.zeros((len(data_set),3),dtype=float)
    for i in range(0,len(data_set)):
        in_file=data_set[i] + "\\analyzed_data\\molting_cycle_data.p"
        tmp=pickle.load( open(in_file, "rb" ) )
        n_molt=1+tmp
        dt_hatch=get_datetime_for_timepoint(data_set[i]+'\\',n_molt[1])
        for j in range(0,3):
            dt_ecd=get_datetime_for_timepoint(data_set[i]+'\\',n_molt[j+2])
            diff=dt_ecd-dt_hatch
            T_ecd[i,j]=(diff.days * 24.) + (diff.seconds/3600.)
    return(T_ecd)

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

# data sets for wild-type (WT) and unc-6
    
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

# to calculate DTC positions from scratch, uncomment this:
               
## get position data
#data=[]
#data.append(get_DTC_position(WT_data_set,correct_for_movement=True))
#data.append(get_DTC_position(unc6_data_set,correct_for_movement=True))
#pickle.dump(data, open('DTC_data_movement_corrected.p', 'wb'))
## get ecdysis data
#ecd_data=[]
#ecd_data.append(get_ecdysis_time_after_L1_ecdysis(WT_data_set))
#ecd_data.append(get_ecdysis_time_after_L1_ecdysis(unc6_data_set))
#pickle.dump(ecd_data, open('ecdysis_data.p', 'wb'))

# and, to calculate DTC positions from scratch, comment this out:
data = pickle.load(open('figure_data\\DTC_data_movement_corrected.p', 'rb'))
ecd_data = pickle.load(open('figure_data\\ecdysis_data.p', 'rb'))

# for converting pixels to micrometers
um_per_pixel=0.1625

# create figure
R=0.45
W=11
fig=plt.figure(2,figsize=(W,R*W))
plt.clf()

q=plt.axes([0,0,1.,1.])
q.axis('off')
q.set_xlim([0,1])
q.set_ylim([0,1])

# create subpanel axes for panel B
W_a=.25
H_a= W_a/10/R
dh_a=0.01
a=[]
for i in range(0,11):
    a.append(plt.axes([.03, .70-(H_a+dh_a)*i, W_a, H_a]))

# create axes for panels C and D
W=0.275    
H=0.35
dh=0.125
b=[]
for i in range(0,2):
    b.append(plt.axes([.38, 0.6-(H+dh)*(i), W,H]))
    
# create axes for panels E and F
W=0.245    
H=0.165
dh=0.03
c=[]
for i in range(0,2):
    c.append(plt.axes([.745, 0.79-(H+dh)*i, W,H]))

# create axes for panels G and H
d=[]
for i in range(0,2):
    d.append(plt.axes([.745, 0.32-(H+dh)*i, W,H]))

### panel B

# load subpanels images
indir=WT_data_set[8]+'\\'
n_molt= get_molting_frame([indir])
n_list = np.array([46,53,61,68,76,84,89,99,107,115,123])
# for each subpanel, get time after L1 ecdysis
T=np.zeros(len(n_list))
for i in range(0,len(n_list)):
    filename=indir + 'analyzed_data\\488nm_straight_%03d.tif' % n_list[i]
    tmp=tif.imread(filename)
    
    if i==0:
        im_data=np.zeros((tmp.shape[0],tmp.shape[1],len(n_list)))    
        datetime0=get_datetime_for_timepoint(indir,n_list[i])
    im_data[:,:,i]=tmp

    datetime=get_datetime_for_timepoint(indir,n_list[i])
    # calculate difference
    diff=datetime-datetime0
    # convert to hours
    T[i]=(diff.days * 24.) + (diff.seconds/3600.)

# plot lines to indicate which images belong to which larval stage
for i in range(1,4):
    ind1=np.where(n_list==n_molt[0][i])[0][0]
    ind2=np.where(n_list==n_molt[0][i+1])[0][0]
    q.plot([.3,.3],[0.70-(H_a+dh_a)*(ind1)+H_a,.70-(H_a+dh_a)*(ind2-1)],'-k')

# plot subpanels
for i in range(0,len(n_list)):
    tmp=im_data[:,:,i]
    a[i].imshow(tmp,cmap='gray',interpolation='nearest',aspect='auto')
    a[i].plot([0,tmp.shape[1]],((tmp.shape[0]-1)/2)*np.ones(2),':w')
    a[i].set_xlim([0,tmp.shape[1]])
    a[i].set_ylim([tmp.shape[0],0])
    a[i].text(-105,65,"%2dh" % np.round(T[i]), color='k', fontsize=10)
    a[i].axis('off')
    if i==0:
        # this is the scale factor in make_straightened_animal_images
        f_sc=2.
        a[i].plot([20,20+20./f_sc/um_per_pixel],[tmp.shape[0]-20,tmp.shape[0]-20],'w',lw=2)

### panel C + D
special_data_set=[[8],[5,2,4]]
special_col='rgb'
cc=0
# plot A-P position versus time
for i in range(0,2):
    # reorder animals so that those that are emphasized by color are plotted
    # last. Start with order in which loaded
    ind_plot=range(0,len(data[i]))
    # remove all animals in special data set
    for j in special_data_set[i]:
        ind_plot.remove(j)
    # and add those to the back
    ind_plot=ind_plot+special_data_set[i]
    for j in ind_plot:
        if (i==0) and (j in special_data_set[i]):
            col_A=[0,0,0]
            col_P=[0,0,0]
            lw=1.5
        elif (i==1) and (j in special_data_set[i]):
            col_A=special_col[cc]
            col_P=special_col[cc]
            cc+=1
            lw=1.5
        else:
            col_A=[0.7,0.7,0.7]
            col_P=[0.7,0.7,0.7]
            lw=1
        b[i].plot(um_per_pixel*data[i][j].s[:,0],data[i][j].T/60,'-', color=col_A,linewidth=lw)
        b[i].plot(um_per_pixel*data[i][j].s[:,1],data[i][j].T/60,'-',color=col_P,linewidth=lw)
        
    # plot histogram of times of ecdysis
    T_ecd=ecd_data[i]
    (h,x)=np.histogram(T_ecd.flatten(),np.arange(0,30,1))
    b[i].barh(x[0:-1],-15*h,left=350,linewidth=0,color='r',height=0.5)

    # set limits and labels
    b[i].set_ylim([30,0])
    b[i].set_xlim([-350,350])
    b[i].plot([0,0],[-1,30],'--k')
    b[i].set_ylabel('Hours after L1 ecdysis',fontsize=12)
    b[i].set_xlabel('Position along A-P axis (mm)')
    b[i].tick_params(axis='both', which='major', labelsize=10)

### panel E+F
# plot D-V position versus time
cc=0
for i in range(0,2):
    # reorder animals so that those that are emphasized by color are plotted
    # last. Start with order in which loaded
    ind_plot=range(0,len(data[i]))
    # remove all animals in special data set
    for j in special_data_set[i]:
        ind_plot.remove(j)
    # and add those to the back
    ind_plot=ind_plot+special_data_set[i]
    for j in ind_plot:
        if (i==0) and (j in special_data_set[i]):
            col_A=[0,0,0]
            col_P=[0,0,0]
            lw=1.5
        elif (i==1) and (j in special_data_set[i]):
            col_A=special_col[cc]
            col_P=special_col[cc]
            cc+=1
            lw=1.5
        else:
            col_A=[0.7,0.7,0.7]
            col_P=[0.7,0.7,0.7]
            lw=1
        c[i].plot(data[i][j].T/60, um_per_pixel*data[i][j].t[:,0], '-', color=col_A,linewidth=lw)
        c[i].plot(data[i][j].T/60, um_per_pixel*data[i][j].t[:,1], '-',color=col_P,linewidth=lw)

    # set limits and labels
    c[i].set_xlim([0,30])
    c[i].set_ylim([-15,20])
    c[i].plot([-1,30],[0,0],'--k')
    plt.setp(c[i], yticks=[-10,0,10,20])
    if i==1:
        c[i].set_ylabel('Position along D-V axis (mm)')
        c[i].set_xlabel('Hours after L1 ecdysis')
    c[i].tick_params(axis='both', which='major', labelsize=10)

## panels G and H - DTC velocity

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
    MT=get_molting_frame(ds)
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
        # check if it has turned, i.e. did not move beyond limit from center of body
        if np.abs(data[i][j].s[-1,k])<(150./um_per_pixel):
            # if not, add to "turn"list
            avg_ind_unc6_turn.append([i,j,k]) 
        else:
            # else, save in "no turn" list...
            if k==0:
                # ... either for DTC.a
                avg_ind_unc6_no_turn_A.append([i,j,k])
            else:
                # ... or else for DTC.p
                avg_ind_unc6_no_turn_P.append([i,j,k])
# collate the index list of the three classes (WT, unc-6 turn and unc-6 no turn DTC.a) in one list
all_ind_list=[avg_ind_WT,avg_ind_unc6_turn,avg_ind_unc6_no_turn_A]

avg_data=[]
# for all three classes calculate the average A-P velocity
for m in range(0,3):
    # calculate A-P velocities
    (t, v) = get_velocity_vs_time(all_ind_list[m],data,DT_L1L3,T_range,DT_acq,W_savg,axis="AP")
    # and get average over all DTCs in class
    (t_avg, v_avg, v_stderr) = get_average_velocity_vs_time(t,v)
    avg_data.append([t_avg,v_avg,v_stderr])    

# calculate dorsal_ventral velocities v_DV for WT
(t, v_DV) = get_velocity_vs_time(all_ind_list[0],data,DT_L1L3,T_range,DT_acq,W_savg,axis="DV")
# get average <v_DV> over all DTCs
(t_avg, v_DV_avg, v_DV_stderr) = get_average_velocity_vs_time(t,v_DV)

### panel G
# plot average A-P and D-V velocity for WT DTCs

# draw lines
d[0].plot([0,0],[-1,1],'--k')
d[0].plot(T_range/60.,[0,0],'--k')

# get WT A-P data
m=0
t_avg=avg_data[m][0]
v_avg=avg_data[m][1]
v_stderr=avg_data[m][2]
# and plot
d[0].plot(t_avg/60.,v_avg,'-',color='r')
d[0].errorbar(t_avg[::3]/60., v_avg[::3], yerr=v_stderr[::3], fmt='.',color='r')

# plot D-V data
d[0].plot(t_avg/60.,2*v_DV_avg,'-',color='k')
d[0].errorbar(t_avg[::3]/60., 2*v_DV_avg[::3], yerr=2*v_DV_stderr[::3], fmt='.',color='k')

# set limits and labels
d[0].set_xlim([-15,10])
d[0].set_ylim([-0.25,0.35])
d[0].set_ylabel('<v> (mm/h)')
plt.setp(d[0], yticks=[-0.2,0,0.2])
plt.setp(d[0], xticks=[-15,-10,-5,0,5,10])
d[0].tick_params(axis='both', which='major', labelsize=10)

### panel H
# plot average A-P velocity for WT DTCs, unc-6 turn DTCs and unc-6 no turn DTC.a
col='rgb'
# plot lines
d[1].plot(T_range/60.,[0,0],'--k')
d[1].plot([0,0],[-1,1],'--k')
# plot A-P veolcities for the three different classes
for m in range(0,3):
    t_avg=avg_data[m][0]
    v_avg=avg_data[m][1]
    v_stderr=avg_data[m][2]
    d[1].plot(t_avg/60.,v_avg,'-',color=col[m])
    d[1].errorbar(t_avg[::3]/60., v_avg[::3], yerr=v_stderr[::3], fmt='.',color=col[m])
# set limits and labels
d[1].set_xlim([-15,10])
d[1].set_ylim([-0.3,0.5])
d[1].set_ylabel('<v_A-P> (um/h)')
plt.setp(d[1], yticks=[-0.2,0,0.2,0.4])
plt.setp(d[1], xticks=[-15,-10,-5,0,5,10])
d[1].tick_params(axis='both', which='major', labelsize=10)
d[1].set_xlabel('Hours after L3 ecdysis')

# print some information about what is plotted. Some of this is needed for the article text

print "Animals plotted with emphasizing colors:"
for i in range(0,2):
    for j in special_data_set[i]:
        if i==0:
            print WT_data_set[j]
        else:
            print unc6_data_set[j]
print ""

# print standard deviation in final DTC position
i=0
N=len(data[i])
s_final=np.zeros((N,2))
for j in range(0,N):
    s_final[j,:]=data[i][j].s[-1,:]
s_std=np.std(um_per_pixel*s_final,axis=0)
lbl='ap'
for i in range(0,2):
    print "DTC.%c, sigma_A-P=%f um" % (lbl[i],s_std[i])
    
# get indeces of unc-6 DTCS, split into those that turn around and those that don't
# for those that don't turn, separate in DTC.a and DTC.p 
i=1
N_tot=0
N_ventral=0
N_no_turn=0
# for each animal
for j in range(0,len(data[i])):
    # for each DTC
    for k in range(0,2):
        N_tot+=1
        if np.abs(data[i][j].s[-1,k])<(150./um_per_pixel):
            N_no_turn += 1
        r=np.where(data[i][j].t[:,k]>0)
        f=float(len(r[0]))/float(len(data[i][j].t[:,k]))
        if f<0.1:
            N_ventral +=1
print "No turn DTCs: %d/%d animals" % (N_no_turn,N_tot)
print "Only ventral DTCs: %d/%d animals" % (N_ventral,N_tot)
