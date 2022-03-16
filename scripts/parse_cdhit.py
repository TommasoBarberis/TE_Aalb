import sys

print(sys.argv)
cdhit_file = "data/TE_copies.fasta.hit.clstr"
output = "test/"

with open(cdhit_file, "r") as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith(">"):
            filename = output+"cluster_"+line.replace(">","").split()[1]+".lst"
            with open(filename, "a") as clust:
                pass
        else:
            seq_name = line.split(">")[1].split("...")[0]
            with open(filename, "a") as clust:
                clust.write(seq_name+"\n")
        if line.startswith(">Cluster 3"):
            break
