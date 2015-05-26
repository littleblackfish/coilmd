These are some simple anlysis scripts. An understanding of the preferred folder hierarchy is necerrary to practically use loopy ones.

## melt.awk

calculates dissociation ratio. 
input is *bubbles.dat*
outputs two floating points, first is the dissociation ratio, second is std error. 
associated temperature can be extracted from current dir name. 
output is usually piped to a file called *melt.dat* which then can be plotted via plot-eb.gp

`melt.awk bubbles.dat >> ../melt.dat`

*loop-melt.sh* is a bash loop that uses this script to generate melt.dat for a full set of temperatures. It should be used from the root directory and takes 2 parameters that are N and type. 

`loop-melt.sh 50 cir`

## size.awk

calculates bubble size distribution.
input is *bubbles.dat*
outputs a non-normalized (count) histogram. 
output is usually piped to a file *bubbledist.dat*

`size.awk bubbles.dat > bubblesize.dat`

## Persistence

**[persistence.py](persistence.py)** is a python script to calculate bond correlations among the polymer axis. It requires [MDAnalysis](http://www.mdanalysis.org/) python package.

It reads a **dcd type trajectory** and a **pdb type topology** file. These can be generated from a vtf type trajectory using **[convertDCD.vmd](convertDCD.vmd)** vmd script. 

`persistence.py topol.pdb traj.dcd`

Generates a file called *persistence.dat*.

**[persistence.gp](persistence.gp)** plots this file and fits an analytical function to extract lambda (pitch) and lp (persistence length). It might be necessary to change the initial conditions in the script to obtain a reasonable fit. 

## plotting scripts

All plotting scripts make use of gnuplot.

* **[plot.gp](plot.gp)** plots all 2 column files matching given pattern using lines.

* **[plot-lp.gp](plot-lp.gp)** plots all 2 column files matching given pattern using lines and points.

* **[plot-eb.gp ](plot-eb.gp)** plots all 3-column files that match the given pattern, using 3rd column for error bars.

* **[plot-kymograph.gp](plot-kymograph.gp)** plots a plain text matrix file (such as bubbles.dat) as a kymograph. Input must be piped in such as : `plot-kymograph.gp < 100/lin/0.18/bubbles.dat`


