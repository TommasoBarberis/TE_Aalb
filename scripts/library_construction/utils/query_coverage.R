#!/usr/bin/Rscript

# https://github.com/clemgoub/TE-Aid/blob/master/consensus2genome.R
# https://github.com/clemgoub/TE-Aid/blob/master/Run-c2g.R

Args 		    =	commandArgs()
blast_file		=	Args[6]
cons_len        =   as.numeric(Args[7])
out_file        =   as.character(Args[8])

cons_coverage=function(blast_file=NULL, cons_len=NULL, out_file=NULL){
    blast=read.table(blast_file, sep="\t")

    #make the coverage matrix
    coverage=matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
    for(i in 1:length(blast$V1)){
        coverage[i,]<-c(rep(0,blast$V6[i]-1),rep(1,abs(blast$V7[i]-blast$V6[i])+1), rep(0,as.numeric(cons_len)-blast$V7[i]))
    }
    coverage<-colSums(coverage)
    write.table(coverage, file="coverage.txt", col.names=F, row.names=F)

}

cons_coverage(blast_file = blast_file,
    cons_len = cons_len,
    out_file = out_file
)