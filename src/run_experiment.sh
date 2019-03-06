#!/bin/bash

for i in `seq 1 10`;
do
    echo Running experiment $i
    GSL_RNG_SEED=$i ./sampler ts.csv > $1/$i.csv

    grep unique $1/$i.csv | cut -d' ' -f8 > $1/mu$i.csv

    grep largest $1/$i.csv | cut -d' ' -f5 > $1/presentations$i.csv
done

