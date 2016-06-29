import numpy as np

class POI:
    def __init__(self,x,lbl):
        self.x=np.array(x)
        self.lbl=lbl

def find_animal_center(POI_list):
    # find slice that contains all three POIs 1,2/v and 3
    Z0=-1
    # first, find all z positions of POIs
    z=np.zeros(len(POI_list),dtype=int)
    for i in range(0,len(POI_list)):
        z[i]=int(POI_list[i].x[2])
    # then, for each z position    
    for i in np.unique(z):
        lbl_list=[]
        # find all labels of POIs with that z position
        for j in range(0,len(POI_list)):
            if int(POI_list[j].x[2])==i:
                lbl_list.append(POI_list[j].lbl)
        # check whether all relevant POIs are in this slice
        if ('1' in lbl_list) and ('3' in lbl_list) and (('2' in lbl_list) or ('v' in lbl_list)):
            # if so, set center slice to this slice
            Z0=i
    return Z0