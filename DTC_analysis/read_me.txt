## Image analysis ##

1) run *<make_small_trans_movie.py> to get a scaled-down trans light movie for annotating ecdyses

2) run *<annotate_molting_cycle.py> to find time of hatching and ecdyses

NOTE: use *<print_molting_cycle_report.py> to get an overview of the timepoints for hatching and ecdyses

3) run <annotate_POI.py> to get POIs 1, 2/v , 3 and d

NOTE: POIs correspond to different anatomical markers: 1 - posterior pharyngeal bulb, 2 - gonad, use this if vulva not yet induced, v - center of vulva, 3 - posterior end of gut and start of anus, d - distal tip cell. There should be only a single slice that contains the POIs 1,2/v and 3. Additional POIs with label 1 ,2/v and 3 can be present in other slices and are used to estimate movement of the animal between slices. 

4) run <annotate_worm_axis.py> to get worm axis

## analysis of DTC dynamics ##

- The function <get_DTC_position.py> can be used to get the position of the DTCs in the (s,t)-coordinate system formed by the animal's A-P and D-V axis.

## article figures ##

- <make_DTC_figure.py> creates Figure 3 in the main text
- <make_DTC_SI_figure.py> creates Supplementary Figure 3 

## Visualization of DTC dynamics ##

- <make_straightened_animal_images.py> is used to create images of straightened animals for each timepoint
- <make_DTC_movie.py> uses images of straightened animals to create movie stills that show the time and central body axis. In addition, the colors are adjusted over the entire duration of DTC migration to correct for increased GFP fluorescence. Still images are assembled into movies using imageJ

## general libraries
- <POI.py> - definition of points of interest (POI) data structure
- <spline_projection.py> - functions to get (s,t) coordinates of cells from (x,y) coordinates and animal's body axis
 
* = check to make final version
