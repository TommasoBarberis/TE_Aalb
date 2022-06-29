"""
Author: Tommaso BARBERIS
Date: 26/05/2022
Description: script allowing to attribute super-family classification
"""

import argparse

def merge_classif(rm_file, rp_file, rc_file):
    
    rm_dict = {}
    with open(rm_file, 'r') as f:
        data = f.readlines()
        for line in data:
            fields = line.split()
            rm_dict[fields[0]] = (fields[3], fields[4])
    
    rp_dict = {}
    with open(rp_file, 'r') as f:
        f.readline()
        data = f.readlines()
        for line in data:
            fields = line.split()
            if "/" in fields[2]:
                classif = fields[2].split("/")
            else:
                classif = [fields[2], fields[3]]
            rp_dict[fields[0]] = (classif[0], classif[1])
    
    rc_dict = {}
    with open(rc_file, 'r') as f:
        data = f.readlines()
        rc_dict = {}
        for line in data:
            fields = line.split()
            if "/" in fields[1]:
                classif = fields[1].split("/")
            else:
                classif = [fields[1], fields[1]]
            rc_dict[fields[0]] = (classif[0], classif[1])

    with open("merged_classif.tsv", "w") as f:
        for seq in rc_dict.keys():
            if seq in rm_dict.keys():
                rm_sup = rm_dict[seq][1]
                rm_inf = rm_dict[seq][0]
            else:
                rm_sup = "-"
                rm_inf = "-"
            if seq in rp_dict.keys():
                rp_sup = rp_dict[seq][1].split("/")[0]
                rp_inf = rp_dict[seq][0]
            else:
                rp_sup = "-"
                rp_inf = "-"
            rc_sup = rc_dict[seq][0]
            rc_inf = rc_dict[seq][1]
            
            # thin tunning
            if rc_sup == "RC":
                rc_sup = "DNA"
            if "Lian" in rp_inf:
                rp_inf = "Lian"
            if "Copia" in rp_inf:
                rp_inf = "Copia"
            if "Lian" in rm_inf:
                rm_sup = "LINE"
            if "Loner" in rm_inf:
                rm_sup = "LINE"
            if any(term in rm_inf for term in ['Bel', 'BEL', 'PAO', 'Pao']):
                rm_inf = "Bel/Pao"
            if any(term in rp_inf for term in ['Bel', 'BEL', 'PAO', 'Pao']):
                rp_inf = "Bel/Pao"
            if any(term in rc_inf for term in ['Bel', 'BEL', 'PAO', 'Pao']):
                rc_inf = "Bel/Pao"


            tmp_line = seq + "\t" + rm_sup + "\t" + rm_inf + "\t" + rp_sup + "\t" + rp_inf + "\t" + rc_sup + "\t" + rc_inf + "\n"
            count = tmp_line.count("-\t-")
            if count < 2:
                sup = "-"
                inf = "-"
                if rm_sup == rp_sup or rm_sup == rc_sup:
                    sup = rm_sup
                elif rp_sup == rm_sup or rp_sup == rc_sup:
                    sup = rp_sup
                elif rc_sup == rm_sup or rc_sup == rp_sup:
                    sup = rc_sup

                
                if "-" in rm_inf and rm_inf != "-":
                    rm_inf = rm_inf.split("-")[0]
                if "_" in rm_inf:
                    rm_inf = rm_inf.split("_")[0]
                if "-" in rp_inf and rp_inf != "-":
                    rp_inf = rp_inf.split("-")[0]
                if "_" in rp_inf:
                    rp_inf = rp_inf.split("_")[0]
                if "-" in rc_inf and rc_inf != "-":
                    rc_inf = rc_inf.split("-")[0]
                if "_" in rm_inf:
                    rc_inf = rc_inf.split("_")[0]

                
                if rm_inf == rp_inf or rm_inf == rc_inf:
                    inf = rm_inf
                if rp_inf == rm_inf or rp_inf == rc_inf:
                    inf = rp_inf
                if rc_inf == rm_inf or rc_inf == rp_inf:
                    inf = rc_inf

                if inf == "-":
                    if rc_inf.startswith(rm_inf) or rm_inf.startswith(rc_inf):
                        inf = rc_inf
                    if rm_inf.startswith(rp_inf) or rp_inf.startswith(rm_inf):
                        inf = rm_inf
                    if rp_inf.startswith(rc_inf) or rc_inf.startswith(rp_inf):
                        inf = rc_inf
                
                if sup == '-':
                    if inf == 'CR1' or inf == 'Dong' or inf == 'I':
                        sup = "LINE"
                    elif  inf == 'Helitron':
                        sup = "DNA"
                    elif inf == 'Penelope':
                        sup = "PLE"
                    elif inf == 'tRNA':
                        sup = "tRNA"
                    
                if sup != '-' or inf != '-':
                    line = (seq + "\t" + sup + "\t" + inf + "\n")
                    f.write(line)

        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-rm", metavar='', help=('''
    filtered output of RepeatMasker with RepBase
    '''))
    parser.add_argument("-rp", metavar='', help=('''
    filtered output of RepeatProteinMask 
    '''))
    parser.add_argument("-rc", metavar='', help=('''
    filtered output of RepeatClassifier
    '''))

    args = parser.parse_args()

    # assign parameters
    rm_file = args.rm
    rp_file = args.rp
    rc_file = args.rc

    if rm_file and rp_file and rc_file:
        merge_classif(rm_file, rp_file, rc_file) 
