#!/usr/bin/gnuplot -p

#set term pdfcairo enhanced size 1024 ,768
#set term pngcairo enhanced
#set out 'bubbles.png'

set key off
set autoscale fix
unset colorbox
#set xlabel 'time'
#set ylabel 'i'
#set cbrange [0:1]

plot '/dev/stdin' matrix w image

