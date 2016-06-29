import numpy as np
import scipy.interpolate as interp
import scipy.optimize as opt

def get_axis_spline(ax):
    # get arc length of axis
    s=np.zeros(ax.shape[0])
    for i in range(1,ax.shape[0]):
        dx=ax[i,:]-ax[i-1,:]
        s[i]=s[i-1]+np.sqrt(np.sum(dx**2))
    # then get spline as function of arc length
    f_spline=1e1
    spline_x = interp.UnivariateSpline(s, ax[:, 0],s=f_spline*ax.shape[0])
    spline_y = interp.UnivariateSpline(s, ax[:, 1],s=f_spline*ax.shape[0])
    # and return it
    return (s[-1],spline_x,spline_y)

def get_spline_unit_vectors(s,spline_x,spline_y):
    # get unit vector along spline, based on derivative at s
    tmp=np.array([spline_x.derivatives(s)[1],spline_y.derivatives(s)[1]])
    e_s=tmp/np.sqrt(np.sum(tmp**2))
    # not get othogonal unit vector
    e_t=np.array([e_s[1], -e_s[0]])
    # CAUTION:: the orientation of e_t with respect to e_s is not uniquely defined!!
    # This depends on the orientation of the gonad (POI 2 or v) with respect to the
    # A-P axis, which need to be determined later!!!!!!
    return (e_s,e_t)

def find_distsq(s,x,spline_x,spline_y):
    # find distance^2 between point on spline with arclength s and point x
    x_s=np.array([spline_x(s), spline_y(s)])
    return np.sum((x_s-x)**2)  

# this is ugly and can be done better, but I am getting tired...
def find_global_min(fun,s,args):
    # get arguments for minimization function <fun>
    x,spline_x,spline_y = args
    # first, find the point that minimizes <fun> in the list of grid point s_i
    min_val=6e66
    min_ind=-1
    # loop over all points s_i
    for i in range(0,len(s)):
        # get function value at point s_i
        val = fun(s[i],x,spline_x,spline_y)
        if val<min_val:
            # if smaller than current smallest value, save value and index i
            min_val=val
            min_ind=i
    # determine range for local minimization, if s_i is the smallest grid point
    # then use range (s_i-1,s_i+1)     
    sb=[-1,-1]
    # deter
    if min_ind==0:
        # unless i=0, then use (0,s_i+1)
        sb[0]=s[0]
    else:
        sb[0]=s[min_ind-1]        
    if min_ind==len(s)-1:
        # or unless i is last grid point, then use (s_i-1,s_i)
        sb[1]=s[-1]
    else:
        sb[1]=s[min_ind+1]        
    # then, determine minimizing arc length s within the range (s_i-1,s_i+1)    
    s_min = opt.fminbound(fun, sb[0], sb[1], args=(x,spline_x,spline_y))
    return s_min

def get_spline_coordinates_st_for_x(x,spline_x,spline_y,L):
    # find spline arclength s that minimzes distance to point x
    # I use a special routine to find the global minimum, if the worm is curved 
    # in a particular way there might be multiple local minima and fminbound
    # can settle on a minimum that is not global!
    s=find_global_min(find_distsq,np.linspace(0,L,100),(x,spline_x,spline_y))
    x_s=np.array([spline_x(s), spline_y(s)])
    # get unit vectors e_s and e_t, along spline and othogonal to spline
    (e_s,e_t) = get_spline_unit_vectors(s,spline_x,spline_y)
    # get the value of t, by projection along e_t 
    # first, calculate the vectorial distance between point at s and position of POI
    dx=x-x_s
    t=np.dot(dx,e_t)
#    # test plot     
#    DX=t*e_t
#    plt.plot([x_s[0],x_s[0]+DX[0]],[x_s[1],x_s[1]+DX[1]], '--y')
    return (s,t)