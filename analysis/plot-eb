#!/bin/bash

# gnuplot scrtipt to plot 3 column files using error bars

command="gnuplot -p -e \"plot "

n=1

for y in $@
do
#  echo $y
  for x in $y*
  do
    command="$command '$x' w errorbars lt $n , '$x' w l lt $n t '',"  
    ((n+=1))
  done

done
command="${command%,}\""

echo $command
#echo $@

eval $command

