library(haven)
library(readxl)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(parallel)
spss_meta <- read_sav("/restricted/projectnb/cteseq/projects/CTE_data_processing/09.28.2023_dataset/CTE RNAseq phenotype data.sav")
csv_meta <- read.csv("/restricted/projectnb/cteseq/projects/CTE_data_processing/all_meta.csv")
counts <- read.csv("/restricted/projectnb/cteseq/projects/CTE_data_processing/all_counts.csv")
AT8_ex_data <- read.csv("/restricted/projectnb/cteseq/jrose/cte_DLFCWSAT8.csv")[ ,2:3]
cte_missing_data <- read_excel("/restricted/projectnb/cteseq/projects/CTE_data_processing/09.28.2023_dataset/CTE CASES FOR MISSING APOE AND TMEM.xlsx")
#Fixes a mismatched sample name
spss_meta$subjid <- gsub("BVAX81", "BVAX081", spss_meta$subjid)
#Changes format to match other rs value
spss_meta$rs3173615 <- gsub("C", "C:C", spss_meta$rs3173615)
spss_meta$rs3173615 <- gsub("G", "G:G", spss_meta$rs3173615)
spss_meta$rs3173615 <- gsub("C:CG:G", "C:G", spss_meta$rs3173615)
#renamed for ease of merging
colnames(cte_missing_data)[1] = "subjid"
colnames(cte_missing_data)[2] = "apoe"
colnames(cte_missing_data)[3] = "rs3173615"
#Merge cte_missing_data into spss_meta file.
spss_merged <- merge(spss_meta, cte_missing_data, by = "subjid", all.x = TRUE)
spss_merged$rs3173615 <- ifelse(spss_merged$rs3173615.x == "" & spss_merged$rs3173615.y != "", spss_merged$rs3173615.y, spss_merged$rs3173615.x)
spss_merged$apoe <- ifelse(is.na(spss_merged$apoe.x) & spss_merged$apoe.y != "", spss_merged$apoe.y, spss_merged$apoe.x)
spss_merged <- subset(spss_merged, select = -c(apoe.x, apoe.y, rs3173615.x, rs3173615.y))
#Merge in additional AT8 data
AT8_ex_data$subjid <- gsub("K", "K-", AT8_ex_data$subjid)
AT8_ex_data$subjid <- gsub("SLI", "SLI-", AT8_ex_data$subjid)
AT8_ex_data$subjid <- gsub("SLI-0", "SLI-", AT8_ex_data$subjid)
spss_merged <- merge(spss_merged, AT8_ex_data, by = "subjid")


#Sets counts gene names (X) to row names, then removing that column
rownames(counts) <- counts$X
counts <- subset(counts, select = -X)
#Filters out gene names that have all zeros
nonzero_genes <- rowSums(counts)!=0
filtered_zeros <- counts[nonzero_genes,]
#Filters out gene names that have 150 or more zeros
ltohf_genes <- rowSums(counts == 0) < 150
filtered_counts <- counts[ltohf_genes,]
#Sorts column names alphabetically to align with 
filtered_counts <- filtered_counts[, order(names(filtered_counts))]



#substitutes K values with only 3 digits after "K-" to have a 0 in front of them
csv_meta$Sample_ID <- gsub("^K-(\\d{3})", "K-0\\1", csv_meta$Sample_ID)
#Does the same as above, but for SampleName instead of Sample_ID.
csv_meta$SampleName <- gsub("^K-(\\d{3}$)", "K-0\\1", csv_meta$SampleName)
#changes Corde_ID in CSV to line up with counts file. Example: ST_AN_K-688 to ST_AN_K.688
csv_meta$Core_ID <- gsub("^ST_AN_K-", "ST_AN_K.", csv_meta$Core_ID)
csv_meta$Core_ID <- gsub("^ST_AN_SLI-", "ST_AN_SLI.", csv_meta$Core_ID)

#Changing the value of 7 to NA in npbraak. It seems odd with a continuous model and only has one person
spss_merged$npbraak <- gsub("7", NA, spss_merged$npbraak)

#Merges meta data files
csv_meta <- csv_meta[,c("SampleName", "RIN", "Core_ID")]
merged_meta <- merge(csv_meta, spss_merged, by.x = "SampleName", by.y = "subjid")

#Transforming apoe, Group/CTE and tmem
merged_meta$DLFC.WS.AT8..cell.mm2 <- gsub("#N/A", NA, merged_meta$DLFC.WS.AT8..cell.mm2)
merged_meta$Group_de <- ifelse(merged_meta$Group %in% c(1,2), "low", ifelse(merged_meta$Group == 3, "high", NA))
merged_meta$apoe_de <- ifelse(merged_meta$apoe %in% c(24, 34, 44), 1, 0)
merged_meta$apoe_de <- ifelse(merged_meta$apoe %in% c(24, 34, 44), 1, 0)
merged_meta$TMEM106B_dom <- ifelse(merged_meta$rs3173615 %in% c("C:C", "C:G"), 1, ifelse(merged_meta$rs3173615 == "G:G", 0, NA))
merged_meta$TMEM106B_invrec <- ifelse(merged_meta$rs3173615 == "C:C", 1, ifelse(merged_meta$rs3173615 %in% c("G:G", "C:G"), 0, NA))



#Sorts Core_ID alphabetically
merged_meta <- merged_meta[order(merged_meta[["Core_ID"]]), ]
#Sets Core_ID to row names for deseq2 #rownames(merged_meta) <- merged_meta$Core_ID

#AT8log
merged_meta$AT8sulcusLog <- log(merged_meta$AT8sulcus)
merged_meta$AT8crestLog <- log(merged_meta$AT8crest)


#Turns AT8 Extra data into numeric instead of characters
merged_meta$DLFC.WS.AT8..cell.mm2 <- as.numeric(merged_meta$DLFC.WS.AT8..cell.mm2)
merged_meta$AT8_total <- merged_meta$DLFC.WS.AT8..cell.mm2
merged_meta <- subset(merged_meta, select = -c(DLFC.WS.AT8..cell.mm2))
#rearrange factor levels. Only binary category to have reversed factors by default. Others are fine
merged_meta$Group_de <- factor(merged_meta$Group_de, levels = c("low", "high"))

#Filtering out NA values on design columns
merged_meta <- subset(merged_meta, !is.na(merged_meta$agedeath))
merged_meta <- subset(merged_meta, !is.na(merged_meta$RIN))
#Filtering counts to only contain Core_ID (columns in counts) in the merged_AT8 df
#counts_filtered <- counts[, names(counts) %in% merged_AT8$Core_ID]
rounded <- round(filtered_counts)



col_continuous <- c("footyrs", "AFE", "totyrs", "agecogsx", "cogdecage", "mincogage", "AT8sulcusLog", "disdur", "DLFC.WS.AT8..cell.mm2", "AT8crestLog")

col_categorical <- c("DementiaHx", "npold", "npoldd", "nppath", "npftdtau", "PathAD", "PathFTD", "CTE", "suicide", "Group_de", "apoe_de", "TMEM106B_dom", "TMEM106B_invrec")

col_ordinal <- c("npavas", "npwmr", "npbraak", "npneur", "npadnc", "npdiff", "npamy", "nparter")

col_categorical_polychotomous <- c("PathLBD", "CTEStage", "Group", "cod", "nplbod")

col_exclude <- c("SampleName", "RIN", "agedeath", "Core_ID", "race", "nphemo", "npold1", "npold2", "npold3", "npold4", "npoldd1", "npoldd2", "npoldd3", "npoldd4", "nppath6", "nppick", "npftdt2", "npcort", "npprog", "npftdt5", "npftdt8", "npftdt9", "npftdt10", "npftdtdp", "micdorfront", "micinfpar", "micalc", "locratio", "subratio", "CA1ratio", "CA23ratio", "CA4ratio", "sport", "rs1990622", "rs3173615", "apoe", "PathPrion", "npoftd", "micsuptemp", "PathMND", "AT8sulcus", "AT8crest")



run_deseq <- mclapply(
  col_categorical,
  function(col_name) {
    function_meta <- merged_meta[!is.na(merged_meta[[col_name]]), ]
    function_meta[[col_name]] <- as.factor(function_meta[[col_name]])
    function_counts <- rounded[, names(rounded) %in% function_meta$Core_ID]
    design_form <- formula(paste("~ agedeath + RIN +", col_name))
    dds <- DESeqDataSetFromMatrix(countData = function_counts,
                                  colData = function_meta,
                                  design = design_form)
    dds <- DESeq(dds)
    return(dds)
  },
  mc.cores = 32
)
names(run_deseq) <- col_categorical
saveRDS(run_deseq, "/restricted/projectnb/cteseq/jrose/cte_deseq_output_cat_binary.rds")
