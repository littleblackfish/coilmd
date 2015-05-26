These are some simple anlysis scripts

# melt.awk

calculates dissociation ratio. 
input is *bubbles.dat*
outputs two floating points, first is the dissociation ratio, second is std error. 
associated temperature can be extracted from current dir name. 
output is usually piped to a file called *melt.dat* which then can be plotted via plot-eb.gp

`melt.awk bubbles.dat >> ../melt.dat`

*loop-melt.sh* uses this script to generate melt.dat for a full set of temperatures.


# size.awk

calculates bubble size distribution.
input is *bubbles.dat*
outputs a non-normalized histogram. 
output is usually piped to a file *bubbledist.dat*

`size.awk bubbles.dat > bubblesize.dat`

