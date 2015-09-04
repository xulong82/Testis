library(xlsx)
library(dplyr)
library(Biobase)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis")
load("data/dataList.rdt")
for(obj in names(dataList)) assign(obj, dataList[[obj]])
source("~/Dropbox/X/function.R")

# Cell type enrichment
Spgonia <- read.delim("meiotic/Spgonia", stringsAsFactors = F)
PL <- read.delim("meiotic/PL", stringsAsFactors = F)
EL <- read.delim("meiotic/EL", stringsAsFactors = F)
LL_Z <- read.delim("meiotic/LL_Z", stringsAsFactors = F)
LL_Z_EP <- read.delim("meiotic/LL_Z_EP", stringsAsFactors = F)
EP <- read.delim("meiotic/EP", stringsAsFactors = F)
EP_LP_D <- read.delim("meiotic/EP_LP_D", stringsAsFactors = F)
LP_D <- read.delim("meiotic/LP_D", stringsAsFactors = F)

cells <- list(Spgonia = Spgonia, EL = EL, LL_Z = LL_Z, LL_Z_EP = LL_Z_EP, EP = EP, EP_LP_D = EP_LP_D, LP_D = LP_D)
sapply(cells, nrow)

myhyper <- function(g1, g2) {  # Hypergeometric
  bg = rownames(gene) # pool in gene
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bg, g2)), length(g1))
} # Pr(count >= length(intersect(g1, g2)))

(input <- rownames(nonp_up))
lapply(cells, function(x) intersect(input, x$gene_name))
sapply(cells, function(x) length(intersect(input, x$gene_name))) / length(input)

number = sapply(cells, function(x) length(intersect(input, x$gene_name))) 
pval = sapply(cells, function(x) myhyper(input, x$gene_name))
pval["EP"] = 1e-12

g = data.frame(cell = names(cells), number = number, pval = pval)

pdf("Pdf/meiotic_cells.pdf", width = 7, height = 4)
ggplot(g, aes(x = cell, y = -log10(pval))) + 
  geom_point(aes(size = number)) + ylim(c(0, 13)) +
  theme_bw() + xlab("") + ylab("-log10(pvalue)")
dev.off()

# EP genes? what's its expression pattern?
ep_genes <- intersect(input, cells$EP$gene_name)
mei_group <- gsub("_.*", "", colnames(meiotic)) 
mei_stage <- unique(mei_group)
mei_group_fac <- factor(mei_group, levels = mei_stage)
meiotic_mean <- sapply(mei_stage, function(x) rowMeans(meiotic[, mei_group == x]))

ep_meiotic <- meiotic_mean[intersect(ep_genes, rownames(meiotic_mean)), ]
ep_meiotic <- ep_meiotic - ep_meiotic[, 1]

pdf("Pdf/meiotic_pattern.pdf", width = 7, height = 5)
plot(1:6, rep(0, 6), type = "l", col = "blue", xlab = "", xaxt = "n", ylab = "", ylim = c(-2, 9))
for(i in 1:nrow(ep_meiotic)) lines(ep_meiotic[i, ], type = "l", col = "grey30")
lines(colMeans(ep_meiotic), type = "b", lwd = 3, col = "red")
axis(1, at=1:6, labels = paste(c(8, 10, 12, 14, 16, 18), "d"))
dev.off()

# EP cells that gradually increase were stalled in translation with Eif4g3 mutation
# What we see in staining EP cells?

# Could it be happening in the CB?
file <- "CB/Meikar_Supp_Tables_revised_XL.xlsx"

cb_protein <- read.xlsx(file, sheetName = "Table S1", startRow = 3, endRow = 91, stringsAsFactors = F) 
cb_protein_name <- lapply(cb_protein$Protein.short.names, function(idx) unlist(strsplit(idx, split = ";"))) 
cb_protein_name <- unique(unlist(cb_protein_name))
cb_protein_name <- gsub("^ ", "", cb_protein_name)

cb_protein_gene <- gsub(",", ";", gsub(" ", "", cb_protein$Gene.Symbol))
cb_protein_gene <- lapply(cb_protein_gene, function(idx) unlist(strsplit(idx, split = ";"))) 
cb_protein_gene <- unique(unlist(cb_protein_gene))

cb_rna <- read.xlsx(file, sheetName = "Table S5", startRow = 3, endRow = 4979, stringsAsFactors = F) 
cb_rna$FPKM <- as.numeric(cb_rna$FPKM)
cb_rna <- cb_rna[cb_rna$FPKM > 30, ]
cb_rna_gene <- cb_rna$symbol

(cb_protein_nonp <- intersect(rownames(nonp_up), cb_protein_gene))
myhyper(rownames(nonp_up), cb_protein_gene)

(cb_rna_nonp <- intersect(rownames(nonp_up), cb_rna_gene))
myhyper(rownames(nonp_up), cb_rna_gene)

write.xlsx(cb_protein_nonp, file = "Chromatoid/intersect.xlsx", sheetName = "CB_Protein_NONP", append = T)
write.xlsx(cb_rna_nonp, file = "Chromatoid/intersect.xlsx", sheetName = "CB_RNA_NONP", append = T)

# data-driving, but not Hspa2-driving analysis
gene <- tpm[biomart$ensembl_transcript_id, ]
gene <- apply(gene, 2, function(x) tapply(x, biomart$external_gene_name, sum))
gene <- gene[! rowMax(gene) > 1e4, ]  # outlier
gene <- sweep(gene, 2, colSums(gene), "/") * 1e6 
gene <- gene[rowMax(gene) > 20, ] 

grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
grp2 <- factor(gsub("[123]", "", colnames(gene)), levels = grp1)
geno <- factor(gsub("[123].*", "", colnames(gene)), levels = c("W", "M"))
frac <- factor(gsub(".*[123]", "", colnames(gene)), levels = c("IN", "PLM", "NONP"))

fit <- apply(log2(gene + 1), 1, function (x) lm(x ~ frac + geno + frac * geno))

summary(fit[["Hspa2"]])
summary(fit[["Hspa8"]])

fit.e = t(sapply(fit, function(x) summary(x)$coefficients[-1, "Estimate"]))
fit.p = t(sapply(fit, function(x) summary(x)$coefficients[-1, "Pr(>|t|)"]))

# Fraction-specific effect of the Eif4g3 mutation

poly = names(fit)[fit.p[, "fracPLM:genoM"] < 0.05]
nonp = names(fit)[fit.p[, "fracNONP:genoM"] < 0.05]
intersect(nonp, names(fit)[fit.e[, "fracNONP:genoM"] > 0])

(mutAll = names(which(rowSums((fit.p < 0.05)[, c("genoM", "fracPLM:genoM", "fracNONP:genoM")]) > 0)))
(nonpAll = unique(c(setdiff(mutAll, poly), nonp)))

sapply(cells, function(x) myhyper(nonp, x$gene_name))
sapply(cells, function(x) myhyper(mutAll, x$gene_name))
sapply(cells, function(x) myhyper(nonpAll, x$gene_name))

hspa2 = names(which(rowSums((fit.pval < 0.05)[, c("fracNONP", "fracNONP:genoM")]) == 2))
hspa2 = intersect(hspa2, names(which(fit.esti[, "fracNONP:genoM"] > 0)))

gk = mmGK(hspa2)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

# sequence motif of the 285 non-poly genes (homer)
