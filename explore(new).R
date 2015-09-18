library(xlsx)
library(dplyr)
library(quantro)
library(Biobase)
library(biomaRt)
library(RColorBrewer)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis")
load("data/dataList.rdt")
for(obj in names(dataList)) assign(obj, dataList[[obj]])
source("~/Dropbox/X/function.R")

myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(geneAll, g2)), length(g1))
} # Pr(count >= length(intersect(g1, g2)))

gene <- tpm[biomart$ensembl_transcript_id, ]
gene <- apply(gene, 2, function(x) tapply(x, biomart$external_gene_name, sum))
gene <- gene[! rowMax(gene) > 1e4, ]  # outlier
gene <- sweep(gene, 2, colSums(gene), "/") * 1e6 
geneAll <- rownames(gene)
gene <- gene[rowMax(gene) > 20, ] 

grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
gene <- gene[, order(factor(gsub("[123]", "", colnames(gene)), levels = grp1))]
grp2 <- factor(gsub("[123]", "", colnames(gene)), levels = grp1)
matboxplot(as.data.frame(log2(gene + 1)), groupFactor = as.numeric(grp2))

geno <- factor(gsub("[123].*", "", colnames(gene)), levels = c("W", "M"))
frac <- factor(gsub(".*[123]", "", colnames(gene)), levels = c("IN", "PLM", "NONP"))
fit <- apply(log2(gene + 1), 1, function (x) lm(x ~ frac + geno + frac * geno))

summary(fit[["Hspa2"]])
summary(fit[["Hspa8"]])

(hspa8 = log2(gene["Hspa8", ] + 1))
(hspa8 = sapply(grp1, function(x) mean(hspa8[grp2 == x])))
(coefs = summary(fit[["Hspa8"]])$coefficients[, "Estimate"])

cols <- brewer.pal(6,"Dark2")
plot(1:6, hspa8, xaxt="n", pch=17, cex=2, xlim=c(0, 7), ylim=c(-1, 12), xlab="")
axis(1, at=1:6, labels=names(hspa2))
abline(h=0)
arrows(1, 0, 1, coefs[1], lwd=3, col=cols[1]) # WIN as baseline
abline(h=coefs[1], col=cols[1])
arrows(2, coefs[1], 2, coefs[1]+coefs[4], lwd=3, col=cols[4]) # MIN - WIN
arrows(3, coefs[1], 3, coefs[1]+coefs[3], lwd=3, col=cols[3]) # WNONP - WIN
segments(3, coefs[1]+coefs[3], 4,coefs[1]+coefs[3], col=cols[3]) # 
arrows(4, coefs[1]+coefs[3], 4, coefs[1]+coefs[3]+coefs[4], lwd=3, col=cols[4])
arrows(4, coefs[1]+coefs[3]+coefs[4], 4, coefs[1]+coefs[3]+coefs[4]+coefs[6], lwd=3, col=cols[6])
arrows(5, coefs[1], 5, coefs[1]+coefs[2], lwd=3, col=cols[2]) # WNONP - WIN
segments(5, coefs[1]+coefs[2], 6,coefs[1]+coefs[2], col=cols[2]) # 
arrows(6, coefs[1]+coefs[2], 6, coefs[1]+coefs[2]+coefs[4], lwd=3, col=cols[4])
arrows(6, coefs[1]+coefs[2]+coefs[4], 6, coefs[1]+coefs[2]+coefs[4]+coefs[5], lwd=3, col=cols[5])
legend("right", names(coefs), fill=cols, bg="white")

fit.e = t(sapply(fit, function(x) summary(x)$coefficients[-1, "Estimate"]))
fit.p = t(sapply(fit, function(x) summary(x)$coefficients[-1, "Pr(>|t|)"]))

# Fraction-specific effect of the Eif4g3 mutation

poly = names(fit)[fit.p[, "fracPLM:genoM"] < 0.05]
nonp = names(fit)[fit.p[, "fracNONP:genoM"] < 0.05]
intersect(nonp, names(fit)[fit.e[, "fracNONP:genoM"] > 0])

(mutAll = names(which(rowSums((fit.p < 0.05)[, c("genoM", "fracPLM:genoM", "fracNONP:genoM")]) > 0)))
(nonpAll = names(which(rowSums((fit.p < 0.05)[, c("genoM", "fracNONP:genoM")]) > 0)))

sapply(cells, function(x) myhyper(nonp, x$gene_name))
sapply(cells, function(x) myhyper(mutAll, x$gene_name))
sapply(cells, function(x) myhyper(nonpAll, x$gene_name))

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("external_gene_name", "chromosome_name"), "external_gene_name", values = nonp, mart)

by_chr = sapply(c(1:19, "X"), function(x) biomart$external_gene_name[biomart$chromosome_name == x])
sapply(by_chr, function(x) paste(x, collapse = ", "))

by_chr_number = sapply(by_chr, length)

pdf("Pdf/by_chr.pdf", width = 9, height = 6)
barplot(by_chr_number)
abline(0, 0)
dev.off()

hspa2 = names(which(rowSums((fit.pval < 0.05)[, c("fracNONP", "fracNONP:genoM")]) == 2))
hspa2 = intersect(hspa2, names(which(fit.esti[, "fracNONP:genoM"] > 0)))

gk = mmGK(hspa2)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

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


# sequence motif of the 285 non-poly genes (homer)
