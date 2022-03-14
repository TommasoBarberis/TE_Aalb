import bamnostic as bs

TE_list = []
bam = bs.AlignmentFile("data/bam/edta_cons_check.sorted.bam", 'rb')
with bam as in_bam:
    for read in in_bam:
        TE_list.append(read.reference_name)