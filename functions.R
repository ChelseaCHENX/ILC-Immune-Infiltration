library(dplyr)
library(tidyr)
library(readxl)
library(matrixStats)
library(grid)

library(ggplot2)
library(ggrepel) 
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# library(GSVA)
library(survminer)
require(survival)
library(gtools)
library(forestplot)
library(data.table)


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


coxConvert = function(x){ x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                        }



                        