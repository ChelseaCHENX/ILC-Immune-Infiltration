rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'

setwd(rootdir)
source('codes/functions.R')

vectra = read.csv("data/heatmap/heatmap_vectra.csv", row.names = 1, check.names = F) # modified per updated meta info

dat0 = read.csv('data/classifier/ilc_classes.csv', row.names = 1, check.names = F)
dat = read.csv("data/heatmap/Meta_Subtype_57genes_BulkRNAseqCohort.csv", row.names = 1, check.names=F) # log2tpm_sel_class2

# changed ILC subtype to be the original 3 TCGA subtypes
dat <- dat %>%
  mutate(`TCGA ILC Subtype` = case_when(
    pred == 0 ~ 'Reactive-like',
    pred == 1 ~ 'Immune-related',
    pred == 2 ~ 'Proliferative',
    TRUE ~ NA  # Handles any unexpected values
  ))
table(dat$`TCGA ILC Subtype`)

dat$PAM50 = factor(dat$PAM50, levels = c('LumA','LumB','Her2','Basal','Normal'))
dat$`mIHC Availability` <- ifelse(rownames(dat) %in% vectra$`Patient ID`, 'Y', 'N')

# gene expression in TPM
gene_cols = c('ADH1B', 'ADH1C', 'ALDH1L1', 'AQP7P1', 'BBOX1', 'C2orf40',
       'CAPN6', 'CD300LG', 'CHRDL1', 'CLDN19', 'CNTFR', 'COL17A1', 'CRYAB',
       'DLK2', 'EFCAB1', 'FABP4', 'FAT2', 'FFAR3', 'FIGF', 'G0S2', 'GFAP',
       'GPX3', 'GRIA4', 'HIF3A', 'ITIH5', 'KCNIP2', 'KLB', 'KLHL13', 'KLK5',
       'KY', 'LGR6', 'LIPE', 'LPL', 'MEOX1', 'MIA', 'NKPD1', 'PAK7', 'PCOLCE2',
       'PDE3B', 'PFKFB1', 'PI16', 'PIK3C2G', 'PLP1', 'RDH5', 'S100B', 'SCARA5',
       'SCN4A', 'SFRP1', 'SLC19A3', 'SLC27A6', 'SOX10', 'STAC2', 'STX11',
       'TMEM132C', 'TNMD', 'TNXB', 'TP63')
tpm = dat[, gene_cols]

# meta data
meta_cols = c("ILC Subtype","TCGA ILC Subtype",'PAM50',"Age at diagnosis","Race","Postmenopausal","Pathological grade","Pathological stage","Ki67","ER","PR","HER2","ER IHC score",'mIHC Availability')
meta = dat[, meta_cols]
## hard coded colname changes (upper / lower case consistency with plotting)
colnames(meta) = c("ILC Subtype", "TCGA ILC Subtype","PAM50","Age at Diagnosis", "Race", "Postmenopausal", "Pathological Grade", "Pathological Stage", "Ki67", "ER", "PR", "HER2", "ER IHC Score", "mIHC Availability")
## change values to factor / numerical
for (col in setdiff(colnames(meta), c("Age at Diagnosis","ER IHC Score","Ki67"))){meta[,col] = as.factor(meta[,col])}
for (col in c("Age at Diagnosis","ER IHC Score","Ki67")){meta[,col]=as.integer(meta[,col])}


# plotting effect settings (complex heatmap)
ftsz = 24
common_legend_params <- list(
  `Age at Diagnosis` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Race = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Postmenopausal = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Grade` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `Pathological Stage` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  Ki67 = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  ER = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ER IHC Score` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PR = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  HER2 = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  PAM50 = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `ILC Subtype` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `TCGA ILC Subtype` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz)),
  `mIHC Availability` = list(title_gp = gpar(fontsize = ftsz, fontface = "bold"), labels_gp = gpar(fontsize = ftsz))
)


colors_er_pr_her2 <- c("Positive" = "#e27ad6", "Negative" = "#8da0cb", 'None'='gray')
col_fun_ki67 = colorRamp2(c(min(meta$Ki67, na.rm = TRUE), max(meta$Ki67, na.rm = TRUE)), c("white", "maroon"))
col_er_score = colorRamp2(c(min(meta$`ER IHC Score`, na.rm = TRUE), max(meta$`ER IHC Score`, na.rm = TRUE)), c("white", "maroon"))
col_age = colorRamp2(c(min(meta$`Age at Diagnosis`, na.rm=TRUE), max(meta$`Age at Diagnosis`, na.rm=TRUE)), c("white", "maroon"))
col_grade = c('1'='pink','2'='magenta','3'='maroon','NA'='gray')
col_race = c('White'='#4878d0', 'Other'='#ee854a')
col_postmenopausal = c('No'='gray', 'Yes'='maroon')
col_stage = c("1A" = "#ADD8E6", "1B" = "#87CEEB", "2A" = "#1E90FF", "2B" = "#0000CD", '3A'='navy', '3C'='black')  
col_ilc_cluster = c('Non-proliferative'='#4878d0', 'Proliferative'='#ee854a')
col_ilc_tcga_cluster = c('Reactive-like'='#66CCCC', 'Immune-related'='#ADD8E6', 'Proliferative'='#ee854a')
col_pam50 = c("LumA" = "#4878d0", "LumB" = "#ee854a", "Her2" = "#6acc64", "Basal" = "#d65f5f", "Normal" = "#956cb4")
col_ihcavail = c('Y'='#ee854a', 'N'='gray')

ha = HeatmapAnnotation(df = meta, simple_anno_size = unit(1, "cm"),
                       annotation_name_gp = gpar(fontsize = 24),
                       annotation_legend_param = common_legend_params,
                       col = list(
                         `TCGA ILC Subtype` = col_ilc_tcga_cluster, 
                         `ILC Subtype` = col_ilc_cluster, 
                         ER = colors_er_pr_her2, PR = colors_er_pr_her2, HER2 = colors_er_pr_her2,
                         Ki67 = col_fun_ki67, `ER IHC Score` = col_er_score, `Age at Diagnosis` = col_age,
                         Race = col_race, Postmenopausal = col_postmenopausal, 
                         `Pathological Stage` = col_stage, `Pathological Grade` = col_grade,
                         PAM50 = col_pam50, `mIHC Availability` = col_ihcavail)
                       )

p = Heatmap(as.data.frame(t(tpm)), top_annotation = ha, 
            heatmap_legend_param = list(title = "log2TPM", at = c(-4, 0, 4), 
                                        title_gp = gpar(fontsize = 24, fontface='bold'), 
                                        labels_gp = gpar(fontsize = 24),
                                        legend_height = unit(5, "cm"),
                                        title_position = "lefttop-rot"
                                        ), 
            height = unit(30, "cm"),column_dend_height = unit(4, "cm"),
            cluster_rows = hclust(dist(t(tpm), method = "euclidean"), method = "ward.D2"), 
            cluster_columns = hclust(dist(tpm, method = "euclidean"), method = "ward.D2"),
            row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 15))


pdf(file.path(rootdir, "images_rebuttal/BulkRNA.pdf"), width=12, height=25)
print(p)
dev.off()