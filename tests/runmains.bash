#!/bin/bash

if [ $# -eq 0 ]; then
    mainlist=`ls ../src/*/*main`
else
    mainlist=`ls $1/src/*/*main`
fi

for x in $mainlist
do
    if [ $# -eq 0 ]; then
	echo "running $x"
    fi
    $x > /dev/null
done
