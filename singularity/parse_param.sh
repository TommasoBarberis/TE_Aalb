#!/bin/bash

## Author: Tommaso Barberis
## Date: 30/03/2022
## Description: script that parse parameters passed to singularity run (or exec) command 
    # -lib: nucleic database for transposable elements, such as RepBase, Dfam or a customized lib

param=($1)
echo $param

## Parsing parameter
for ((index=0; index <= ${#param[@]}; index++)); do
    case ${p} in
    -lib)
        NUCL_BANK="${param[index+1]}"
        shift
        shift
        ;;
    -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
    esac
done

echo $NUCL_BANK