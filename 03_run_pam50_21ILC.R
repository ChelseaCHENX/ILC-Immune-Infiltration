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

# to exclude normal
prob = PAM50Preds$subtype.proba %>% as.data.frame()# %>% select(-Normal)

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

table(mg$PAM50)
table(mg[,c("PAM50","class_2c")])

# write.csv(mg, "data/classifier/2_class_pam50.csv", quote=F)
# Decide not to merge RNA-seq with as pam50 original data, as GSE10886 is microarray

#---------------- HEATMAP

col_ilc_cluster = c('Non-proliferative'='#4878d0', 'Proliferative'='#ee854a')
col_pam50 = c("LumA" = "#4878d0", "LumB" = "#ee854a", "Her2" = "#6acc64", "Basal" = "#d65f5f", "Normal" = "#956cb4")

ha = HeatmapAnnotation(df = mg[,c('PAM50','ILC Subtype')],
                       col = list(`ILC Subtype`=col_ilc_cluster,PAM50=col_pam50))

p = Heatmap(t(as.matrix(mg[,c('LumA','LumB','Her2','Basal','Normal-like')])),
            heatmap_legend_param = list(title = "Probability", 
                                        title_gp = gpar(fontsize = 10, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 10),
                                        legend_height = unit(3, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            top_annotation=ha, height = unit(3, "cm"),width=unit(8,'cm'))


pdf(file.path(wkdir, "images/PAM50_prob.pdf"), width=9, height=4)
print(p)
dev.off()

