"""
Author: Tommaso Barberis

Date: 26/04/2022

Description: parse .distmat file from distmat emboss's tool. It return the average of all distances.

Dependancies:
    -

Tips:
"""

import argparse


def compute_avg(dist_file):
    sum = 0
    nb_dist = 0

    with open(dist_file, "r") as mat:
        line = mat.readline()
        
        c = 0
        while line:
            if c > 7:
                val_lst = line.split("\t")[c-6:-2]
                for val in val_lst:
                    sum += float(val)
                    nb_dist += 1
            
            line = mat.readline()
            c += 1

    print(dist_file + ":\t" + str(sum/nb_dist))


if __name__ == "__main__":    
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", metavar='', help=('''
    cd-hit-est output file.
    '''))
    
    args = parser.parse_args()

    # assign parameters
    distmat_file = args.m 

    compute_avg(distmat_file)
