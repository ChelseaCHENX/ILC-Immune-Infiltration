rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'
setwd(rootdir)
source('codes/functions.R')

# supp table---------------------------------------------------
vectra = read.csv("data/heatmap/heatmap_vectra.csv", row.names = 1, check.names = F) # modified per updated meta info
dat = read.csv("data/heatmap/Meta_Subtype_57genes_BulkRNAseqCohort.csv", row.names = 1, check.names=F) # log2tpm_sel_class2

dat$PAM50 = factor(dat$PAM50, levels = c('LumA','LumB','Her2','Basal','Normal'))
dat$`mIHC Availability` <- ifelse(rownames(dat) %in% vectra$`Patient ID`, 'Yes', 'No')
dat$`RNA-seq Availability` = 'Yes'

dat$`Patient ID` = rownames(dat)
ihc_counter = as.data.frame(table(vectra$`Patient ID`)); colnames(ihc_counter) = c('Patient ID','ROI Number')

mg = merge(dat, ihc_counter, by='Patient ID', all.x=T)
meta_cols = c('Patient ID','Sample ID',"Tumor type","Age at diagnosis","Race","Postmenopausal","Pathological grade","Pathological stage","Ki67","ER","PR","HER2","ER IHC score","PAM50","ILC Subtype",'RNA-seq Availability','mIHC Availability','ROI Number')
meta = mg[,meta_cols]

write.csv(meta, 'data_rebuttal/supp/sample_info.csv', row.names=F)

# Immune pattern annotations---------------------------------------------------
dat1 = read.csv("data/corr/data_nmflabeled.csv", check.names=F, row.names='Sample Name') # log2tpm_sel_class2, has immune infiltration patterns
dat2 = read.csv("data/heatmap/Meta_Subtype_57genes_BulkRNAseqCohort.csv", row.names = 1, check.names=F) # log2tpm_sel_class2, has pam50 and ILC subtypes
dat2$`Patient ID` <- rownames(dat2)

mg = merge(dat1, dat2[,c('Patient ID','PAM50','ILC Subtype')], by.x='Patient ID', by.y='Patient ID', all.x=TRUE)
rownames(mg) <- rownames(dat1)
cols_sorted = c("Patient ID","ILC Subtype","PAM50","Age at diagnosis","Race","Postmenopausal","Pathological grade","Pathological stage","Ki67","ER","PR","HER2","ER IHC score", 'Immune Infiltration Pattern',
'Stromal B cells','Stromal CD4 T cells','Stromal CD8 T cells','Stromal Treg cells','Stromal Macrophages','Tumor B cells','Tumor CD4 T cells','Tumor CD8 T cells','Tumor Treg cells','Tumor Macrophages','Tumor PanCK cells')
mg_sorted = mg[cols_sorted]
## adding supp2, the raw cell density per ROI with meta data
write.csv(mg_sorted, 'data_rebuttal/supp/cell_density_roi.csv', row.names=T)

meta_cols = c("ILC Subtype","PAM50","Age at diagnosis","Race","Postmenopausal","Pathological grade","Pathological stage","Ki67","ER","PR","HER2","ER IHC score", 'Immune Infiltration Pattern')
meta = mg[, meta_cols]

for (col in setdiff(meta_cols, c("Age at diagnosis","ER IHC score","Ki67"))){meta[,col] = as.factor(meta[,col])}
for (col in c("Age at diagnosis","ER IHC score","Ki67")){meta[,col]=as.integer(meta[,col])}

# hard coded
colnames(meta) = c("ILC Subtype","PAM50","Age at Diagnosis", "Race", "Postmenopausal", 
"Pathological Grade", "Pathological Stage", 
"Ki67", "ER", "PR", "HER2", "ER IHC Score", 'Immune Infiltration Pattern')


ftsz = 24
common_legend_params <- list(
  # Histology = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Age at Diagnosis` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Race = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Postmenopausal=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Grade`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Stage`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Ki67=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  ER=list(title_gp = gpar(fontsize = ftsz, fonXtface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ER IHC Score`= list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PR=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  HER2=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PAM50=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ILC Subtype`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),

#   `mIHC Availability`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz))
`Immune Infiltration Pattern`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz))
)

colors_er_pr_her2 <- c("1" = "#e27ad6", "0" = "#8da0cb", 'None'='gray')
col_fun_ki67 = colorRamp2(c(min(meta$Ki67, na.rm = TRUE), max(meta$Ki67, na.rm = TRUE)), c("white", "maroon"))
col_er_score = colorRamp2(c(min(meta$`ER IHC Score`, na.rm = TRUE), max(meta$`ER IHC Score`, na.rm = TRUE)), c("white", "maroon"))
col_age = colorRamp2(c(min(meta$`Age at Diagnosis`, na.rm=TRUE), max(meta$`Age at Diagnosis`, na.rm=TRUE)), c("white", "maroon"))

col_grade = c('1'='pink','2'='magenta','3'='maroon','NA'='gray')
#c('pink', 'maroon','gray'); names(col_grade) = levels(factor(meta$`Pathological Grade`, exclude = NULL))

col_race = c('White'='#4878d0', 'Other'='#ee854a')
col_postmenopausal = c('0'='gray', '1'='maroon')
col_stage = c("1A" = "#ADD8E6", "1B" = "#87CEEB", "2A" = "#1E90FF", "2B" = "#0000CD", '3A'='navy', '3C'='black')  
col_ilc_cluster = c('Non-proliferative'='#4878d0', 'Proliferative'='#ee854a')
col_pam50 = c("LumA" = "#4878d0", "LumB" = "#ee854a", "Her2" = "#6acc64", "Basal" = "#d65f5f", "Normal" = "#956cb4")

# col_ihcavail = c('Y'='#ee854a', 'N'='gray')
set2_colors = c('#A1C9F4','#FFB482','#8DE5A1','#FF9F9B','#D0BBFF')#brewer.pal(5, "Pastel1")
col_immunepattern = col_immunepattern <- c('1' = set2_colors[1], '2' = set2_colors[2], '3' = set2_colors[3], '4' = set2_colors[4], '5' = set2_colors[5])

immune_infiltration_order <- order(meta$`Immune Infiltration Pattern`)
df_ha = meta[immune_infiltration_order,]

ha = HeatmapAnnotation(df = df_ha, simple_anno_size = unit(1, "cm"),
                       annotation_name_gp = gpar(fontsize = 24),
                       annotation_legend_param = common_legend_params,
                       col = list(ER=colors_er_pr_her2, PR=colors_er_pr_her2, HER2=colors_er_pr_her2,
                                  Ki67=col_fun_ki67, `ER IHC Score`=col_er_score, `Age at Diagnosis`=col_age,
                                  Race=col_race,  Postmenopausal=col_postmenopausal, 
                                  `Pathological Stage`=col_stage, `Pathological Grade`=col_grade,
                                  `Immune Infiltration Pattern`=col_immunepattern,
                                  PAM50=col_pam50, `ILC Subtype`=col_ilc_cluster
                                  )
                       )

ordered_columns <- rownames(df_ha)
heatmap_matrix <- matrix(runif(0), nrow = 1, ncol = 97)
colnames(heatmap_matrix) <- ordered_columns  # Ensure columns are named and ordered properly

ht_opt$ROW_ANNO_PADDING = unit(20, "cm")

p = Heatmap(heatmap_matrix, top_annotation = ha, 
            heatmap_legend_param = list(
                                        title_gp = gpar(fontsize = 24, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 24),
                                        legend_height = unit(5, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            height = unit(0, "cm"),
            cluster_columns = FALSE, cluster_rows = FALSE,
            row_names_gp = gpar(fontsize = 24), column_names_gp = gpar(fontsize = 15),
            show_heatmap_legend = FALSE 
            )

pdf(file.path(rootdir, "images_rebuttal/immunePattern_info.pdf"), width=20, height=15)
draw(p, show_annotation_legend = FALSE)
# draw(ha)
# print(p)
dev.off()


# ILC vs immune patterns contingency
dat = read.csv('data_rebuttal/supp/cell_density_roi.csv', check.names=F)
contingency_table <- table(dat$`ILC Subtype`, dat$`Immune Infiltration Pattern`)
print(contingency_table)
chi_square_test <- chisq.test(contingency_table)
print(chi_square_test)

# number of ROIs per sample
table(dat[,c('Patient ID')])
