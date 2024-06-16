library(haven)
library(readxl)
library(dplyr)
library(ggplot2)
library(mice)
library(corrplot)
library(olsrr)
set.seed(123)
#Multi-col: Group, 

merged_meta <- readRDS("/restricted/projectnb/cteseq/jrose/meta_imputation/merged_meta.rds")

chiif_meta <- subset(merged_meta, select = c("SampleName", "chii_f", "sport", "npavas", "nphipscl", "npwmr", "ParknismHx", 
                                             "Group_de", "CTEStage", "npadnc", "nppath", "npbraak", "PathAD", "AT8_total", 
                                             "apoe_de", "apoe", "PathLBD", "agedeath", "npftdtau", "nparter", "npamy", "npdiff",
                                             "npneur", "DementiaHx", "totyrs", "AFE", "Year" ))
rownames(chiif_meta) <- chiif_meta$SampleName
data_chiif <- chiif_meta[,2:ncol(chiif_meta)]

chiif_model <- lm(chii_f ~ ., data = data_chiif)

ols_step_best <- ols_step_best_subset(chiif_model)
saveRDS(ols_step_best, "/restricted/projectnb/cteseq/jrose/meta_imputation/ols_step_best.rds")
