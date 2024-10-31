rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'

setwd(rootdir)
# source('codes/functions.R')
library(ComplexHeatmap)
library(circlize)

vectra = read.csv("data/heatmap/heatmap_vectra.csv", row.names = 1, check.names = F) # modified per updated meta info

mat_cols = c('Stromal B cells','Stromal CD4 T cells','Stromal CD8 T cells','Stromal Treg cells','Stromal Macrophages','Tumor B cells','Tumor CD4 T cells','Tumor CD8 T cells','Tumor Treg cells','Tumor Macrophages')#,'Tumor PanCK cells')
sel_cols = c("ILC Subtype",'PAM50',"Race", "Postmenopausal", "Pathological Grade", "Pathological Stage", "ER", "PR", "HER2")

mat = vectra[,mat_cols]

# meta = vectra[,c(sel_cols,'ER IHC Score','Ki67','Age at Diagnosis')]
meta_cols = c("ILC Subtype","PAM50","Age at diagnosis","Race","Postmenopausal","Pathological Grade","Pathological Stage","Ki67","ER","PR","HER2","ER IHC Score")
meta = vectra[,meta_cols]
# hard coding, editing upper/lower cases
colnames(meta) = c("ILC Subtype","PAM50","Age at Diagnosis","Race","Postmenopausal","Pathological Grade","Pathological Stage","Ki67","ER","PR","HER2","ER IHC Score")


for (col in setdiff(colnames(meta), c("Age at Diagnosis","ER IHC Score","Ki67"))){meta[,col] = as.factor(meta[,col])}
for (col in c("Age at Diagnosis","ER IHC Score","Ki67")){meta[,col]=as.integer(meta[,col])}

# for (col in sel_cols){meta[,col] = as.factor(meta[,col])}
# for (col in c("ER IHC Score",'Ki67')){meta[,col]=as.numeric(meta[,col])}
# for (col in c('Age at Diagnosis')){meta[,col]=as.integer(meta[,col])}
# meta$PAM50 = factor(meta$PAM50, levels = c('LumA','LumB','Her2','Basal','Normal'))

# table(meta$`Age at Diagnosis`)
# table(meta$`ER IHC Score`)
# table(meta$`Ki67`)
# # table(meta$`PAM50`)
# table(meta$`HER2`)
# summary(meta$`Age at Diagnosis`)

ftsz = 24
common_legend_params <- list(
  # Histology = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Age at Diagnosis` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Race = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Postmenopausal=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Grade`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Stage`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Ki67=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  ER=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ER IHC Score`= list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PR=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  HER2=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PAM50=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ILC Subtype`=list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz))
)

colors_er_pr_her2 <- c("Positive" = "#e27ad6", "Negative" = "#0a4eed")
col_fun_ki67 = colorRamp2(c(min(meta$Ki67, na.rm = TRUE), max(meta$Ki67, na.rm = TRUE)), c("white", "maroon"))
col_er_score = colorRamp2(c(min(meta$`ER IHC Score`, na.rm = TRUE), max(meta$`ER IHC Score`, na.rm = TRUE)), c("white", "maroon"))
col_age = colorRamp2(c(min(meta$`Age at Diagnosis`, na.rm=TRUE), max(meta$`Age at Diagnosis`, na.rm=TRUE)), c("white", "maroon"))

col_grade = c('1'='pink','2'='maroon','NA'='gray')
#c('pink', 'maroon','gray'); names(col_grade) = levels(factor(meta$`Pathological Grade`, exclude = NULL))

col_race = c('White'='#4878d0', 'Other'='#ee854a')
col_postmenopausal = c('No'='blue', 'Yes'='maroon')
col_stage = c("1A" = "#ADD8E6", "1B" = "#87CEEB", "2A" = "#1E90FF", "2B" = "#0000CD", '3A'='navy')  
col_ilc_cluster = c('Non-proliferative'='#4878d0', 'Proliferative'='#ee854a')
col_pam50 = c("LumA" = "#4878d0", "LumB" = "#ee854a", "Her2" = "#6acc64", "Basal" = "#d65f5f", "Normal" = "#956cb4")

col_pt = c(
  "TP17_M510" = "#4878d0", "TP17_M802" = "#ee854a", "TP17_M882" = "#6acc64", 
  "TP18_M121" = "#d65f5f", "TP18_M202" = "#956cb4", "TP18_M251" = "#8c613c",
  "TP18_M278" = "#dc7ec0", "TP18_M302" = "#797979", "TP18_M329" = "#d5bb67",
  "TP18_M355" = "#82c6e2", "TP18_M372" = "#b47cc7", "TP18_M403" = "#c4e17f",
  "TP18_M95"  = "#292421"
)

# both
ha1 = HeatmapAnnotation(df = meta, simple_anno_size = unit(.8, "cm"),
                       annotation_name_gp = gpar(fontsize = 20),
                       annotation_legend_param = common_legend_params,
                       col = list(ER=colors_er_pr_her2, PR=colors_er_pr_her2, HER2=colors_er_pr_her2,
                                  Ki67=col_fun_ki67, `ER IHC Score`=col_er_score, `Age at Diagnosis`=col_age,
                                  Race=col_race, Postmenopausal=col_postmenopausal,  
                                  `Pathological Stage`=col_stage, `Pathological Grade`=col_grade,
                                  `ILC Subtype`=col_ilc_cluster, PAM50=col_pam50)
                       )

ha2 = HeatmapAnnotation(Patient = vectra$`Patient ID`, simple_anno_size = unit(.5, "cm"),
                       col = list(Patient=col_pt),
                       annotation_name_gp = gpar(fontsize = 20)
                       )

mat_log <- log2(mat + 1)
rownames(mat_log) = NULL

rows.cor <- cor(mat_log, use = "pairwise.complete.obs", method = "spearman")
hclust.row <- hclust(as.dist(1-rows.cor), method='ward.D2')

p = Heatmap(t(mat_log), top_annotation = ha1, bottom_annotation = ha2,
            heatmap_legend_param = list(title = "log2 Cell Density (/mm2)", 
                                        title_gp = gpar(fontsize = 20, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 20),
                                        legend_height = unit(3, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            height = unit(10, "cm"),width=unit(30,'cm'),column_dend_height = unit(4, "cm"),
            cluster_rows = hclust.row, 
            cluster_columns = hclust(dist(mat_log, method = "euclidean"), method = "ward.D2"),
            row_names_gp = gpar(fontsize = 20), column_names_gp = gpar(fontsize = 20))


pdf(file.path(rootdir, "images_rebuttal/Vectra_Both.pdf"), width=30, height=12)
print(p)
dev.off()

#---------------------------------------------------
# non-proliferative
meta1 = subset(meta, `ILC Subtype`=='Non-proliferative')
vectra1 = subset(vectra, `ILC Subtype`=='Non-proliferative')
mat1= vectra1[,mat_cols]

ha1 = HeatmapAnnotation(df = meta1, simple_anno_size = unit(.8, "cm"),
                       annotation_name_gp = gpar(fontsize = 24),
                       annotation_legend_param = common_legend_params,
                       col = list(ER=colors_er_pr_her2, PR=colors_er_pr_her2, HER2=colors_er_pr_her2,
                                  Ki67=col_fun_ki67, `ER IHC Score`=col_er_score, `Age at Diagnosis`=col_age,
                                  Race=col_race, Postmenopausal=col_postmenopausal,  `ILC Subtype`=col_ilc_cluster,
                                  `Pathological Stage`=col_stage, `Pathological Grade`=col_grade),
                        show_legend=F
                       )

ha2 = HeatmapAnnotation(Patient = vectra1$`Patient ID`, simple_anno_size = unit(.5, "cm"),
                       col = list(Patient=col_pt),
                       annotation_name_gp = gpar(fontsize = 24),
                       show_legend=F
                       )

mat_log <- log2(mat1 + 1); rownames(mat_log) = NULL
p = Heatmap(t(mat_log), top_annotation = ha1, bottom_annotation = ha2,
            heatmap_legend_param = list(title = "log2 Cell Density (/mm2)", at = c(-4, 0, 4), 
                                        title_gp = gpar(fontsize = 18, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 18),
                                        legend_height = unit(3, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            show_heatmap_legend = F,
            height = unit(10, "cm"),width=unit(12,'cm'),column_dend_height = unit(4, "cm"),
            cluster_rows = F,#hclust(dist(t(mat_log), method = "euclidean"), method = "ward.D2"), 
            cluster_columns = hclust(dist(mat_log, method = "euclidean"), method = "ward.D2"),
            row_names_gp = gpar(fontsize = 24), column_names_gp = gpar(fontsize = 10))

pdf(file.path(rootdir, "images/Vectra_nonProlif.pdf"), width=18, height=12)
print(p)
dev.off()

# proliferative
meta1 = subset(meta, `ILC Subtype`=='Proliferative')
vectra1 = subset(vectra, `ILC Subtype`=='Proliferative')
mat1= vectra1[,mat_cols]

ha1 = HeatmapAnnotation(df = meta1, simple_anno_size = unit(.8, "cm"),
                       annotation_name_gp = gpar(fontsize = 24),
                       annotation_legend_param = common_legend_params,
                       col = list(ER=colors_er_pr_her2, PR=colors_er_pr_her2, HER2=colors_er_pr_her2,
                                  Ki67=col_fun_ki67, `ER IHC Score`=col_er_score, `Age at Diagnosis`=col_age,
                                  Race=col_race, Postmenopausal=col_postmenopausal,  `ILC Subtype`=col_ilc_cluster,
                                  `Pathological Stage`=col_stage, `Pathological Grade`=col_grade),
                        show_legend=F
                       )

ha2 = HeatmapAnnotation(Patient = vectra1$`Patient ID`, simple_anno_size = unit(.5, "cm"),
                       col = list(Patient=col_pt),
                       annotation_name_gp = gpar(fontsize = 24),
                       show_legend=F
                       )

mat_log <- log2(mat1 + 1); rownames(mat_log) = NULL
p = Heatmap(t(mat_log), top_annotation = ha1, bottom_annotation = ha2,
            heatmap_legend_param = list(title = "log2 Cell Density (/mm2)", at = c(-4, 0, 4), 
                                        title_gp = gpar(fontsize = 18, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 18),
                                        legend_height = unit(3, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            show_heatmap_legend = F,
            height = unit(10, "cm"),width=unit(12,'cm'),column_dend_height = unit(4, "cm"),
            cluster_rows = F,#hclust(dist(t(mat_log), method = "euclidean"), method = "ward.D2"), 
            cluster_columns = hclust(dist(mat_log, method = "euclidean"), method = "ward.D2"),
            row_names_gp = gpar(fontsize = 24), column_names_gp = gpar(fontsize = 10))


pdf(file.path(rootdir, "images/Vectra_Prolif.pdf"), width=18, height=12)
print(p)
dev.off()


