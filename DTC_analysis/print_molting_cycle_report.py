import numpy as np
import os
import cPickle as pickle

data_dir='W:\\Jeroen\\timelapse_data\\28_08_2015_qIs56\\'

event=['Hatch','L1','L2', 'L3', 'L4']

# find all animals in data directory
tmp=os.listdir(data_dir)
for s in tmp:
    # for each animal
    in_dir=data_dir + s + "\\analyzed_data\\"
    if os.path.isdir(in_dir):
        print s +":",
        in_file=in_dir + "molting_cycle_data.p"
        if os.path.isfile(in_file):
            T=pickle.load( open(in_file, "rb" ) )
            s=""
            for i in range(0,len(T)):
                if T[i] != -1:
                    s=s + event[i] + ":%d " % (T[i]+1)
            if len(s)>0:
                print "\t"+ s
        else:
            print "\tMolting cycle not yet annotated!"