"""
Author: Tommaso BARBERIS
Date: 05/05/2022
Description: parse blastn output, check if hits come from the same family, generate new fasta file.
Previous commands:
    blastn -query consensi.fa -db db/RepBase27.03.multifasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -out consensi_rnd2.o
    awk '{OFS="\t"; if ($3 >= 80 && ($4/$13 > 0.8) && $4 >= 80) {print $1,$2,$3,$4,$4/$13}}' consensi_rnd2.o > consensi_rnd2.filtered.o
"""

import argparse
from Bio import SeqIO

def parse_fasta(fasta):
    """
    Given a fasta file, it returns a dict where the key is the ID of the sequence and the value is the sequence itself.
    It can be used with the fasta file used for cd-hit-est.
    """

    seqs = {} # key: name of the copy, value: sequence
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name = seq_record.description
        seq = seq_record.seq
        seqs[name] = str(seq).upper() # it allows to avoid errors while running Refiner
    return seqs

def parse_tsv(tsv_file, fasta_file):
    """
    Parse tsv file from the filtered blastn output.
    """

    with open(tsv_file, "r") as tsv:
        data = tsv.readlines()
        consensi = [(line.split()[0], line.split()[1]) for line in data]
        hit_lst = list(set(consensi))

        def classified_fasta(lst, fasta_file):
            seq_dict = parse_fasta(fasta_file)

            for cons in lst:
                seq = seq_dict[cons[0]]
                with open(fasta_file + ".classified", "a+") as f:
                    f.write(">" + cons[1] + "\n" + seq + "\n")
                with open("classified", "a+") as f:
                    f.write(cons[0] + "\t" + cons[1] + "\n")


        if len(hit_lst) == len(data): # there aren't ambigous hits
            classified_fasta(hit_lst, fasta_file)

        else: # process ambigous hits
            uniq_lst = []
            uniq_consensi =  {}
            for cons in consensi: # dict[cons name] = [list of db hits]
                if cons[0] not in uniq_consensi.keys():
                    uniq_consensi[cons[0]] = cons[1]
                else:
                    if isinstance(uniq_consensi[cons[0]], list):
                        tmp_list = uniq_consensi[cons[0]]
                    else:
                        tmp_list = [uniq_consensi[cons[0]]]                
                    tmp_list.append(cons[1])
                    uniq_consensi[cons[0]] = tmp_list
    
            for k in uniq_consensi.keys():
                
                if isinstance(uniq_consensi[k], list): # if the are severals hits 

                    def which_sep(string):
                        """
                        Detect the first fields separator that is used in the classification name ('-' or '_')
                        """
                        
                        for c in string:                        
                            if c == "-":
                                return "-"
                            if c == "_":
                                return "_"
                    
                    # take the prefix of the first member as pattern when possible, else take the full name
                    try:
                        sep = which_sep(uniq_consensi[k][0])                        
                        prefix = uniq_consensi[k][0].split(sep)[0]
                    except:
                        prefix = uniq_consensi[k][0]
                    
                    decision = True
                    # if all members have the same prefix, all hits come from the same family
                    for member in uniq_consensi[k][1:]:                    
                        if member.startswith(prefix[:-1]):
                            pass
                        else:
                            decision = False
                    
                    if decision:
                        uniq_lst.append((k, prefix))
                else:
                    uniq_lst.append((k, uniq_consensi[k]))
            
            # write fasta
            classified_fasta(uniq_lst, fasta_file)
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", metavar='', help=('''
    filtered blastn results file (consensi_rnd2.filtered.o)
    '''))
    parser.add_argument("-f", metavar='', help=('''
    FASTA file with unclassified consensi
    '''))

    args = parser.parse_args()

    # assign parameters
    tsv_file = args.b
    fasta_file = args.f
    
    parse_tsv(tsv_file, fasta_file)