rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'
setwd(rootdir)
source('codes/functions.R')

# METABRIC ER+ --------------------------------------------------------------
## preprocessing
library(data.table)
df = fread('/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/yap/data/normcts/metabric.nolan.expr.csv') %>% as.data.frame()
rownames(df) = df$V1; df$V1=NULL

gex = df

meta = read.csv('../yap/data/metabric/brca_metabric_clinical_data.csv', row.names=1)
pts = intersect(colnames(df), rownames(meta))
meta = meta[pts,]
summary(meta$Age.at.Diagnosis, useNA="ifany")

meta$agedx = as.numeric(meta[pts,]$Age.at.Diagnosis)
meta$stage = as.numeric(meta[pts,]$Tumor.Stage)

gsva = as.data.frame(t(read.csv('data/gsva/metabric_log2cpm_Macrophages.gmt.ssgsea.csv', check.names = F, row.names = 1)))
meta$Sig = as.numeric(gsva[pts,]$ma_neg_pred)

dat = meta
dat[['DSS_time']] = dat[['Relapse.Free.Status..Months.']]
dat[['DSS_status']] = dat[['Relapse.Free.Status']]
dat[['OS_time']] = dat[['Overall.Survival..Months.']]
dat[['OS_status']] = ifelse(dat[['Overall.Survival.Status']] == '0:LIVING', 0, 1)

dat[['ER_status']] = dat[['ER.Status']]
dat[['HER2_status']] = dat[['HER2.Status']]
dat[['PR_status']] = dat[['PR.Status']]
dat[['Tumor_size']] = dat[['Tumor.Size']]
dat[['Tumor_stage']] = dat[['Tumor.Stage']]
dat[['Histology']] = ifelse(dat[['Tumor.Other.Histologic.Subtype']] == 'Lobular', 'ILC', 'other')
dat[['Histology']] = factor(dat[['Histology']], levels=c('other','ILC'))

## DSS, all histology
vars_sel = c('DSS_time', 'DSS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive'
         ) %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()

dim(plot.dat) # 1096

x.cut= surv_cutpoint(plot.dat, time="DSS_time", event="DSS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(DSS_time, DSS_status) ~ Sig+Histology+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) 
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_DFS.csv', quote=FALSE)

km_fit <- survfit(Surv(DSS_time, DSS_status) ~ Sig, data = plot.dat)
p = ggsurvplot(km_fit, data = plot.dat,
           pval = TRUE,            # Show p-value from log-rank test
           conf.int = F,        # Show confidence intervals
           risk.table = F,      # Show the risk table
           ggtheme = theme_minimal() + 
                     theme(plot.title = element_text(hjust = 0.5)), # Center the title
           xlab = "Time (Days)",  # X-axis label
           ylab = "Relapse Free Survival Probability",  # Y-axis label
           title = "METABRIC ER+ Tumors",  # Title of the plot
           legend.title = "TAM-Low Signature",
           palette = c("#6497B1", "#F7766D"))  # Custom legend title
pdf(file.path(rootdir, "images_rebuttal/KM METABRIC DFS.pdf"), width=4, height=3)
print(p)
dev.off()

## DSS, IDC
vars_sel = c('DSS_time', 'DSS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive' &
         `Cancer.Type.Detailed` == 'Breast Invasive Ductal Carcinoma') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()

dim(plot.dat) # 821

x.cut= surv_cutpoint(plot.dat, time="DSS_time", event="DSS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(DSS_time, DSS_status) ~ Sig+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) # no Histology
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_DFS_IDC.csv', quote=FALSE)

## DSS, ILC
vars_sel = c('DSS_time', 'DSS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive' &
         `Cancer.Type.Detailed` == 'Breast Invasive Lobular Carcinoma') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()

dim(plot.dat) # 88

x.cut= surv_cutpoint(plot.dat, time="DSS_time", event="DSS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(DSS_time, DSS_status) ~ Sig+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) # no Histology
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_DFS_ILC.csv', quote=FALSE)


## OS, all histology
vars_sel = c('OS_time', 'OS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive'
         ) %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()
dim(plot.dat) # 1096

x.cut= surv_cutpoint(plot.dat, time="OS_time", event="OS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_time, OS_status) ~ Sig+Histology+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) 
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_OS.csv', quote=FALSE)

km_fit <- survfit(Surv(OS_time, OS_status) ~ Sig, data = plot.dat)
p = ggsurvplot(km_fit, data = plot.dat,
           pval = TRUE,            # Show p-value from log-rank test
           conf.int = F,        # Show confidence intervals
           risk.table = F,      # Show the risk table
           ggtheme = theme_minimal() + 
                     theme(plot.title = element_text(hjust = 0.5)), # Center the title           xlab = "Time in Days",  # X-axis label
           ylab = "Overall Survival Probability",  # Y-axis label
           xlab = 'Time (Days)',
           title = "METABRIC ER+ Tumors",  # Title of the plot
           legend.title = "TAM-Low Signature",
           palette = c("#6497B1", "#F7766D"))  # Custom legend title
pdf(file.path(rootdir, "images_rebuttal/KM METABRIC OS.pdf"), width=4, height=3)
print(p)
dev.off()

## OS, IDC
vars_sel = c('OS_time', 'OS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive' &
         `Cancer.Type.Detailed` == 'Breast Invasive Ductal Carcinoma') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()

x.cut= surv_cutpoint(plot.dat, time="OS_time", event="OS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_time, OS_status) ~ Sig+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) # no Histology
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_OS_IDC.csv', quote=FALSE)

## OS, ILC
vars_sel = c('OS_time', 'OS_status','Sig', 'HER2_status','PR_status', 'Histology',  'agedx','Tumor_size','Tumor_stage')

plot.dat <- dat %>%
  filter(`ER.Status` == 'Positive' &
         `Cancer.Type.Detailed` == 'Breast Invasive Lobular Carcinoma') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()
dim(plot.dat)

x.cut= surv_cutpoint(plot.dat, time="OS_time", event="OS_status", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_time, OS_status) ~ Sig+HER2_status+PR_status+agedx+Tumor_size+Tumor_stage, data = plot.dat); print(summary(cox)) # no Histology
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/METABRIC_OS_ILC.csv', quote=FALSE)


# SCANB ER+ --------------------------------------------------------------
## preprocessing
meta = read.csv('/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/yap/data/scanb/scanb_samples_processed.tsv', 
sep='\t', check.names=F, row.names=1)

gsva = as.data.frame(t(read.csv('data/gsva/scanb_log2cpm_Macrophages.gmt.ssgsea.csv', check.names = F, row.names = 1)))

dat = meta[rownames(gsva),]
dat$Sig = gsva$ma_neg_pred



dat[['Histology']] = ifelse(dat[['Histological_Type']] == 'Lobular', 'ILC', 'other')
dat[['Histology']] = factor(dat[['Histology']], levels=c('other','ILC'))

dat[['HER2_status']] = ifelse(dat[['HER2']] == 'POS', 'pos', 'other')
dat[['HER2_status']] = factor(dat[['HER2_status']], levels=c('other','pos'))

dat[['PR_status']] = ifelse(dat[['PgR_1perc']] == 'POS', 'pos', 'other')
dat[['PR_status']] = factor(dat[['PR_status']], levels=c('other','pos'))

dat$NHG <- factor(dat$NHG, levels = c("G1", "G2", "G3"), ordered = TRUE);dat$NHG_numeric <- as.integer(dat$NHG)
dat$NStage <- factor(dat$NStage, levels = c("N0","N1", "N2", "N3"), ordered = TRUE);dat$NStage_numeric <- as.integer(dat$NStage)
dat$TStage <- factor(dat$TStage, levels = c("T0","T1", "T2", "T3", "T4"), ordered = TRUE);dat$TStage_numeric <- as.integer(dat$TStage)

## OS, all histology
vars_sel = c('OS_days', 'OS_event','Sig','HER2_status', 'PR_status','Histology', 'Age','TumorSize','NStage_numeric','MStage')
for (var in vars_sel[3:10]){print(var);print(summary(dat[[var]], useNA="always"))}

plot.dat <- dat %>%
  filter(`ER_1perc` == 'POS'
         ) %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()

dim(plot.dat)

x.cut= surv_cutpoint(plot.dat, time="OS_days", event="OS_event", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_days, OS_event) ~ Sig+Histology+ Age + HER2_status + PR_status + TumorSize + NStage_numeric+MStage, data = plot.dat); print(summary(cox)) 
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/SCANB_OS.csv', quote=FALSE)

km_fit <- survfit(Surv(OS_days, OS_event) ~ Sig, data = plot.dat)
p = ggsurvplot(km_fit, data = plot.dat,
           pval = TRUE,            # Show p-value from log-rank test
           conf.int = F,           # Show confidence intervals
           risk.table = F,         # Show the risk table
           ggtheme = theme_minimal() + 
                     theme(plot.title = element_text(hjust = 0.5)), # Center the title
           xlab = "Time (Days)",  # X-axis label
           ylab = "Overall Survival Probability",  # Y-axis label
           title = "SCAN-B ER+ tumors",  # Title of the plot
           legend.title = "TAM-Low Signature",
           palette = c("#6497B1", "#F7766D"))  # Custom legend title
pdf(file.path(rootdir, "images_rebuttal/KM SCANB OS.pdf"), width=4, height=3)
print(p)
dev.off()

## OS, IDC
vars_sel = c('OS_days', 'OS_event','Sig','HER2_status', 'PR_status','Histology', 'Age','TumorSize','NStage_numeric','MStage')
for (var in vars_sel[3:10]){print(var);print(summary(dat[[var]], useNA="always"))}

plot.dat <- dat %>%
  filter(`ER_1perc` == 'POS' &
         `Histological_Type` == 'Ductal') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()
dim(plot.dat) # 2199

x.cut= surv_cutpoint(plot.dat, time="OS_days", event="OS_event", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_days, OS_event) ~ Sig+ Age + HER2_status + PR_status + TumorSize + NStage_numeric+MStage, data = plot.dat); print(summary(cox)) # No histology
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/SCANB_OS_IDC.csv', quote=FALSE)

## OS, ILC
vars_sel = c('OS_days', 'OS_event','Sig','HER2_status', 'PR_status','Histology', 'Age','TumorSize','NStage_numeric','MStage')

plot.dat <- dat %>%
  filter(`ER_1perc` == 'POS' &
         `Histological_Type` == 'Lobular') %>%
  dplyr::select(all_of(vars_sel)) %>%
  na.omit()
dim(plot.dat) # 373

for (var in vars_sel[3:10]){print(var);print(summary(plot.dat[[var]], useNA="always"))}


x.cut= surv_cutpoint(plot.dat, time="OS_days", event="OS_event", variables="Sig")
x.cat <- surv_categorize(x.cut)
stopifnot(all(rownames(x.cat)==rownames(plot.dat)))
plot.dat$Sig = factor(x.cat$Sig, levels=c('low','high'))

cox <- coxph(Surv(OS_days, OS_event) ~ Sig+ Age + HER2_status + PR_status + TumorSize + NStage_numeric, data = plot.dat); print(summary(cox)) # No histology or MStage (all M0)
res = as.data.frame(summary(cox)$conf.int); res$pval = as.vector(summary(cox)$coefficients[,'Pr(>|z|)'])   
write.csv(res,'data/survival/SCANB_OS_ILC.csv', quote=FALSE)
