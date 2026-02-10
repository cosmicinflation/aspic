#!/bin/bash

mainlist=`ls $1/src/*/*main`
njobs=`nproc`
count=0
for x in $mainlist; do
    echo
    echo '---------------------->'
    echo 'Starting test program: '$x
    echo '---------------------->'
    echo
    ($x) &
    count=$[$count + 1]
    [[ $((count%njobs)) -eq 0 ]] && wait
done
wait
