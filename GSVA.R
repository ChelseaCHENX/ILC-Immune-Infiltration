#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
exprpath = args[1] # exprpath = '../data/gex/tpm.csv'
dbpath = args[2] # dbpath = '../data/geneset/c2.cp.v2023.1.Hs.symbols.gmt'
resPrefix = args[3] # resPrefix = '../data/gsva/ILC21pt_MsigdbCanonical'
method = args[4] # method = 'ssgsea'

print(args)

library(GSVA)

parse_gmt = function(gmt_path){
        db = readLines(gmt_path)
        geneset = list()
        for (line in db){
            line = 
            words = as.vector(strsplit(line, "\\s{1,}")[[1]])
            set_name = words[1]
            genes = words[-c(1,2)]
            geneset[[set_name]] = genes
        }
        return(geneset)}

# GSVA input (RNAseq): gene X cell, normalized counts

# gsva(expr, gset.idx.list, annotation,
#          method=c("gsva", "ssgsea", "zscore", "plage"),
#          kcdf=c("Gaussian", "Poisson", "none"),
#          abs.ranking=FALSE,
#          min.sz=1,
#          max.sz=Inf,
#          parallel.sz=0,
#          parallel.type="SOCK",
#          mx.diff=TRUE,
#          tau=switch(method, gsva=1, ssgsea=0.25, NA),
#          ssgsea.norm=TRUE,
#          verbose=TRUE)

# expr: gene x cell, current version packageVersion('GSVA') >> ‘1.32.0’, does not specifiy 'is rnaseq'

geneset = parse_gmt(dbpath) 
mat = as.matrix(read.csv(exprpath, header=T, row.names=1)) 

# res = GSVA::gsva(mat, geneset, kcdf='Gaussian', method=method, parallel.sz=8) # it is updated since 2024, thus the old one does not work, per 188
ssgsea_par <- ssgseaParam(mat, geneset)  # all other values are default values
res <- gsva(ssgsea_par)

write.csv(res, file=paste(resPrefix,'.',method,'.csv', sep=''), quote=F)
