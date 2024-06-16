---
title: "Seurat tutorial"
output: 
---

First block of code is creting a seurat object, filtering out specific samples
```{r}
library(dplyr)
library(Seurat)
library(patchwork)

#loading datasets
counts_data <- read.csv("/restricted/projectnb/cteseq/projects/ALS_rnaseq/data/raw/exp_counts.csv", row.names = 1)
meta_data <- read.csv("/restricted/projectnb/cteseq/projects/ALS_rnaseq/data/raw/exp_meta.csv")
#filter metadata under Primary_Dx for ALS only 
meta_data_filtered <- filter(meta_data, Primary_Dx == "ALS")
meta_data_filtered <- as.data.frame(meta_data_filtered)
rownames(meta_data_filtered) <- meta_data_filtered$Core_ID
#Gets each unique value from the Core_ID column in metadata and assigns it to "sample"
sample <- unique(meta_data_filtered$Core_ID)
#subsets the values from "sample" from the counts data so only ALS data remains
counts_data_filtered <- counts_data[,sample]
#passing counts and meta data to a seurat objext
obj <- CreateSeuratObject(counts = counts_data_filtered, meta.data = meta_data_filtered)
obj
```


```{r}

```




Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
