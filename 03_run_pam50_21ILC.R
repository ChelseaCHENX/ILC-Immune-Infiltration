library(dplyr)
library(GEOquery)
library(genefu)
data(pam50.robust)

wkdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune/'
setwd(wkdir)
source('codes/functions.R')

gex = read.csv('data/gex/tpm.csv', row.names=1)
targets_genes = data.frame("Gene.Symbol"=rownames(gex))

PAM50Preds<-molecular.subtyping(sbt.model = "pam50", 
data=as.matrix(t(gex)),
annot=targets_genes,do.mapping=F)

prob = PAM50Preds$subtype.proba %>% as.data.frame()

p = Heatmap(prob)

top_prob_column <- apply(prob, 1, function(row) names(which.max(row))) %>% as.data.frame()
colnames(top_prob_column) = c("PAM50")

ilcclass = read.csv("data/classifier/ilc_classes.csv", row.names=1)

stopifnot(all(rownames(top_prob_column)==rownames(ilcclass)))
mg = cbind(top_prob_column, ilcclass$class_2c)
stopifnot(all(rownames(mg)==rownames(prob)))
mg = cbind(mg, prob)
colnames(mg) = c('PAM50','ILC Subtype', 'Basal','Her2','LumA','LumB','Normal-like')
mg = mg[,c('PAM50','ILC Subtype', 'LumA','LumB','Her2','Basal','Normal-like')]


write.csv(mg, "data/classifier/2_class_pam50.csv", quote=F)
