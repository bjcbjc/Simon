#!/usr/bin/bash

fpattern=$1

for fn in $(ls $fpattern); do
    echo $fn
    msg=`samtools view -H $fn | grep EOF`
    #if [ "$msg" != $'\n' ]; then
#	echo $fn
#	echo $msg
#    fi
done
