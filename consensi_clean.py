from Bio import SeqIO

seq_list = []

for seq_record in SeqIO.parse("data/consensi/RM2_consensi.fa", "fasta"):
    seq_list.append((seq_record.id, seq_record.name, seq_record.description, len(seq_record.seq), seq_record.seq))

with open("data/consensi/RM2_duplicate_headers.txt", "r") as f:
    headers = f.readlines()
    
    seq_info = {}
    for header in headers:
        info = header.split()
        seq_id = info[0].replace(">","")
        seq_num = info[1]
        seq_info[seq_id] = seq_num

for seq in seq_list:
    if seq[0] in seq_info.keys():
        new_id = seq[0]+"."+str(seq_info[seq[0]])

        # update the dict
        if seq_info[seq[0]] == 1:
            del seq_info[seq[0]]
        else:
            new_dup = int(seq_info[seq[0]])-1
            seq_info[seq[0]] = new_dup

    else:
        new_id = seq[0]

    # add the modified header to the new file with its sequence
    with open("RM2_consensi_clean.fa", "a") as new_cons:
        new_cons.write(">"+new_id+"\n"+str(seq[4])+"\n")

    

        