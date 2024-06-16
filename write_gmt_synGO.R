library(readxl)


synGO <- read_excel("/restricted/projectnb/cteseq/jrose/SynGO_bulk_download_release_20231201/syngo_ontologies.xlsx")
gene_ids <- synGO$hgnc_symbol
gene_ids <- trimws(gsub(", ", "\t", gene_ids))
GO_names <- synGO$name
GO_names <- trimws(gsub("\\(GO:\\d+\\)", "", GO_names))
GO_ID <- synGO$id

to_gmt <- paste(GO_names, GO_ID, gene_ids, sep = "\t")
writeLines(to_gmt, "synGO.gmt")
