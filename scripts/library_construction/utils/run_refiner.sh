#!/bin/bash

# Author: Tommaso Barberis
# Date: 08/04/2022
# Description: run Refiner on samples

work_dir=$(pwd)
for dir in $(ls -1 -d cluster_*); do 
    cd $work_dir

    path=$(realpath $dir)
    cd $path
    for clst in $(ls -1); do 
        Refiner $clst
    done

done