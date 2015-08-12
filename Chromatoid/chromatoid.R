library(xlsx)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")

load("data/dataList.rdt")
for(obj in names(dataList)) assign(obj, dataList[[obj]])

file <- "Chromatoid/Meikar_Supp_Tables_revised_XL.xlsx"

cb_protein <- read.xlsx(file, sheetName = "Table S1", startRow = 3, endRow = 91, stringsAsFactors = F) 
cb_protein_symbol <- gsub(",", ";", gsub(" ", "", cb_protein$Gene.Symbol))
cb_protein_symbol <- lapply(x, function(idx) unlist(strsplit(idx, split = ";"))) 
cb_protein_symbol <- unique(unlist(cb_protein_symbol))

cb_rna <- read.xlsx(file, sheetName = "Table S5", startRow = 3, endRow = 4979, stringsAsFactors = F) 
cb_rna$FPKM <- as.numeric(cb_rna$FPKM)
cb_rna_symbol <- cb_rna$symbol

(cb_protein_nonp <- intersect(rownames(nonp), cb_protein_symbol))
(cb_rna_nonp <- intersect(rownames(nonp), cb_rna_symbol))

write.xlsx(cb_protein_nonp, file = "Chromatoid/intersect.xlsx", sheetName = "CB_Protein_NONP", append = T)
write.xlsx(cb_rna_nonp, file = "Chromatoid/intersect.xlsx", sheetName = "CB_RNA_NONP", append = T)
