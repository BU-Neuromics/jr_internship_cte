library(parallel)
library(DESeq2)

merged_meta <- readRDS("/restricted/projectnb/cteseq/jrose/meta_imputation/merged_meta.rds")

#!!!categorizing data columns by type!!!
col_continuous <- c("aestot", "footyrs", "AFE", "bistot", "tbri", "tmi", "tgec", "maxaggsum", 
                    "CDStot", "chii_f", "chii_r", "chii_g", "faqtot", "GDStot", "subratio", "CA1ratio", "CA23ratio", "CA4ratio", 
                    "agecogsx", "cogdecage", "mincogage", "disdur", "AT8_total")

col_binary <- c("DementiaHx", "ParknismHx", "nphemo", "npold", "nppath", 
                "nptdpb", "nptdpc", "nptdpe", "npftdtau", "PathAD", "PathFTD", 
                "CTE", "csparCTE", "sleepact", "Group_de", "TMEM106B_dom", "TMEM106B_invrec")
for (i in col_binary) {
  merged_meta[[i]] <- as.factor(merged_meta[[i]])
}

col_ordinal <- c("CTEStage", "npavas", "npwmr", "npbraak", "npneur", "npadnc", "npdiff", "npamy", "nparter")

for (i in col_ordinal) {
  merged_meta[[i]] <- as.factor(merged_meta[[i]])
}

col_exclude <- c("SampleName", "RIN", "Core_ID", "race", "nphemo", "npold1", "npold2", "npold3", "npold4",
                 "npoldd1", "npoldd2", "npoldd3", "npoldd4", "nppath6", "nppick", "npftdt2", "npcort", 
                 "npprog", "npftdt5", "npftdt8", "npftdt9", "npftdt10", "npftdtdp", "micdorfront", "micinfpar", 
                 "micalc", "locratio", "subratio", "CA1ratio", "CA23ratio", "CA4ratio", "sport", "rs1990622", 
                 "rs3173615", "apoe", "PathPrion", "npoftd", "micsuptemp", "PathMND", "AT8sulcus", "AT8crest", "chii_r")

#categorical <- c(apoe, race, nphipscl, nplbod, PathLBD, PathMND, sport, Group, rs1990622, rs3173615)

#not_using <- c(nppath6, nppick, npftdt2, npcort, npprog, npftdt5, npftdt8, npftdt9, npftdt10, npftdtdp, npoftd, PathPrion, sleepcb, suicide, Batch, Year, npoldd, npold1, npold2, npold3, npold4, npoldd1, npoldd2, npoldd3, npoldd4, locratio)
col_combined <- c(col_binary, col_continuous, col_ordinal)


run_deseq <- mclapply(
  col_combined,
  function(i) {
    print(i)
    function_meta <- merged_meta[!is.na(merged_meta[[i]]), ]
    function_counts <- rounded[, names(rounded) %in% function_meta$Core_ID]
    design_form <- formula(paste("~ agedeath + RIN + Year +", i))
    dds <- DESeqDataSetFromMatrix(countData = function_counts,
                                  colData = function_meta,
                                  design = design_form)
    dds <- DESeq(dds)
    return(dds)
  },
  mc.cores = 16
)
names(run_deseq) <- col_combined
saveRDS(run_deseq, "/restricted/projectnb/cteseq/jrose/cte_corr_work/cte_deseq_output_for_corr.rds")

