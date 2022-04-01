"""
Author: Tommaso Barberis
Date: 01/04/2022
Description: custom script to obtain hierarchy file for PopoolationTE2
Input:file with a header by line for each consensus
"""

import sys

with open('te-hierarchy.txt', "a+") as h:
    h.write("id\tfamily\torder\n")

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    for line in lines:
        fields = line.split("#")
        te_id = fields[0]
        fields = fields[1].split(" ")[0]
        if "/" in fields:
            fields = fields.split("/")
            order = fields[0]
            family = fields[1]
        else:
            order = fields
            family = fields
        
        
        with open('te-hierarchy.txt', "a+") as h:
            h.write(te_id + '\t' + family + '\t' + order + "\n")        