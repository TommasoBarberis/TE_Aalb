#!/usr/bin/Rscript

# https://github.com/clemgoub/TE-Aid/blob/master/consensus2genome.R
# https://github.com/clemgoub/TE-Aid/blob/master/Run-c2g.R

Args 		    =	commandArgs()
blast_file		=	Args[6]
cons_len        =   as.numeric(Args[7])
out_file        =   as.character(Args[8])

cons_coverage=function(blast_file=NULL, cons_len=NULL, out_file=NULL){
    blast = read.table(blast_file, sep="\t")
    new_rownames <- c(paste(blast$V2,"-",as.character(blast$V8),":",as.character(blast$V9), sep = ""))

    #make the coverage matrix
    coverage = matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
    for(i in 1:length(blast$V1)){
        coverage[i,]<-c(rep(0, blast$V6[i]-1),rep(1, abs(blast$V7[i]-blast$V6[i])+1), rep(0, as.numeric(cons_len)-blast$V7[i]))
    }
    
    rownames(coverage) <- new_rownames
    coverage<-colSums(coverage)
    write.table(coverage, file=out_file, col.names=F, row.names=F)

}

cons_coverage(blast_file = blast_file,
    cons_len = cons_len,
    out_file = out_file
)