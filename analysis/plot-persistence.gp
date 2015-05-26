#!/usr/bin/gnuplot

# gnuplot script to plot and fit persistent length data
# lambda is pitch
# lp is persistence length 

set fit quiet
set xlabel '{/Symbol t} (bases)'
set ylabel 'C ({/Symbol t})'


f(x)=exp(-x/lp)*(a+(1-a)*cos(2*pi*x/lambda))
lp=100
lambda=11
a=1


fit f(x) 'persistence' via a,lp,lambda

plot 'persistence' notitle lc 1  pt 7 ps 0.5 , f(x) lc 1 lt 1  lw 4 t 'CG-dsDNA',\

show var l
