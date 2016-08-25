## Image analysis ##

1) run *<downsizeTimelapseImages.py> to get a scaled-down trans light movie for annotating ecdyses

2) create a txt file 'skin.txt' to annotate the ecdyses timepoints

3) run *<GUI_seamCells.py> to manually select the body axis, straighten the images and manually label seam cells

## plot results ##

4) use *<seamCell_results_lineage.py> to seam cell lineage for a single worm (figure 2 panel B)

5) use *<seamCell_results_singleSeamCellMovie.py> to create a movie of a single seam cell (figure 2 panel C)

6) use *<seamCell_results_divTime> to create figure 2 panel D

7) use *<seamCell_results_mistakeMatrix> to create figure 2 panel G

8) use *<seamCell_results_divSequenceMatrix> to create supplementary figure 2

Note: timeLapseFun.py is a general library that contains many useful functions.
