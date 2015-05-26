#!/bin/bash

# script to plot all two column files that match a pattern

command="gnuplot -p -e \"plot "

for y in $@
do
#  echo $y
  for x in $y*
  do
    command="$command '$x' w l,"  
  done

done
command=${command%,}\"

echo $command
#echo $@

eval "$command &"

