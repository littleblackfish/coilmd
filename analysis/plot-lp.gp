#!/bin/bash

# gnuplot script to plot 2 column files using lines AND points

command="gnuplot -p -e \"plot "

for y in $@
do
#  echo $y
  for x in $y*
  do
    command="$command '$x' w lp,"  
  done

done
command=${command%,}\"

echo $command
#echo $@

eval $command

