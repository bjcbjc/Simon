#!/bin/bash

for f in $(ls /nethome/bjchen/Projects/Simon/RNA/varcall/*.vcf); do
    echo $f
    echo `cat $f | grep -v ^# | cut -f9 | sort | uniq | wc -l`
done
