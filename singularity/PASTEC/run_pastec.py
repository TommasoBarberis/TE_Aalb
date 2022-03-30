# -*-coding:Utf-8 -*

"""
Author: Tommaso Barberis

Date: 30/03/2022

Description: wrapper for PASTEC inside the singularity container. It allows to parse some parameters.
"""

import subprocess, time
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-l", "--library", metavar='', help=('''
Path to Repbase database  directory specifically formatted for REPET, containing repbase20.05_ntSeq_cleaned_TE.fa and
repbase20.05_aaSeq_cleaned_TE.fa. The last one is RepBase20.05_REPET.embl.tar.gz.
'''))
parser.add_argument("-i", "--fasta", metavar='', help=('''
input fasta file name [compulsory] [format: fasta]
'''))

args = parser.parse_args()
repbase_dir = args.lib

if __name__ == "__main__":
    start = time.time()
    
    PASTEC_PATH = "/opt/PASTEC_linux-x64-2.0/commons/tools/PASTEClassifier.py"

    command = ['python2.7', PASTEC_PATH] 
    subprocess.run(command)
            
    end = time.time()
    elapsed = end - start
    print(f"\n\n\nElapsed time: {elapsed}")