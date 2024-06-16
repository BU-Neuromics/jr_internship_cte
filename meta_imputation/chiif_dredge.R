library(haven)
library(readxl)
library(dplyr)
library(ggplot2)
library(mice)
library(corrplot)
library(MuMIn)
library(olsrr)
library(glmnet)
library(MASS)
library(car)

set.seed(123)
#global option needs to be set in dredge, essentially requiring complete case data
options(na.action = "na.fail")
#Multi col: Group, CTE, nplbod

merged_meta <- readRDS("/restricted/projectnb/cteseq/jrose/meta_imputation/merged_meta.rds")
#chiif_meta <- subset(merged_meta, select = c("SampleName", "chii_f", "sport", "npavas", "nphipscl", "npwmr", "ParknismHx", "Group_de", "CTEStage", "npadnc", "nppath", "npbraak", "PathAD", "AT8_total", "apoe_de", "apoe", "Group", "CTE", "PathLBD", "agedeath", "npftdtau", "nplbod", "nparter", "npold", "npold1", "npold2", "npold3", "npold4", "npamy", "npdiff", "npneur", "DementiaHx", "totyrs", "AFE", "Year", "Mammillary body", "Insula", "Nucleus Accumbens", "Caduate", "Putamen", "Substantia Innominata", "Optic nerve", "Inferior Frontal Cortex", "Globus Pallidus", "Thalamus", "Cerebellum", "Superior Frontal Cortex", "Motor Cortex", "Calcarine Cortex", "SN", "Anterior Temporal", "Inferior Parital Cortex", "Superior Temporal", "Ent cortex", "Amygdala", "LC", "Hip CA1", "Hip CA2/3", "Hip CA4", "Dorsolateral Frontal Cortex"))
# , "Mammillary body", "Insula", "Nucleus Accumbens", "Caduate", "Putamen", "Substantia Innominata", "Optic nerve", "Inferior Frontal Cortex", "Globus Pallidus", "Thalamus", "Cerebellum", "Superior Frontal Cortex", "Motor Cortex", "Calcarine Cortex", "SN", "Anterior Temporal", "Inferior Parital Cortex", "Superior Temporal", "Ent cortex", "Amygdala", "LC", "Hip CA1", "Hip CA2/3", "Hip CA4", "Dorsolateral Frontal Cortex"))
chiif_meta <- subset(merged_meta, select = c("SampleName", "chii_f", "sport", "npavas", "nphipscl", "npwmr", "ParknismHx", 
                                             "Group_de", "CTEStage", "npadnc", "nppath", "npbraak", "PathAD", "AT8_total", 
                                             "apoe_de", "apoe", "PathLBD", "agedeath", "npftdtau", "nparter", "npamy", "npdiff",
                                             "npneur", "DementiaHx", "totyrs", "AFE", "Year"))


#complete cases required for dredge
cc_meta_chiif <- na.omit(chiif_meta)

rownames(cc_meta_chiif) <- cc_meta_chiif$SampleName
data_chiif <- cc_meta_chiif[,2:ncol(cc_meta_chiif)]

chiif_model <- lm(chii_f ~ ., data = data_chiif)

chiif_dredge <- dredge(chiif_model)
saveRDS(chiif_dredge, "/restricted/projectnb/cteseq/jrose/meta_imputation/chiif_dredge.rds")

