library(dplyr)
library(biomaRt)
library(Biobase)
library(xlsx)
library(genefilter)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")
load("data/ge_tpm.rdt")
tpm <- ge_tpm
# tpm <- ge_tpm[rowMax(as.matrix(ge_tpm)) > 20, ] 
sample <- colnames(tpm)
grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
grp2 <- factor(gsub("[123]", "", sample), levels = grp1)

tpm_nonp <- tpm %>% dplyr::select(contains("NONP"))  # Non-poly
tpm_nonp <- tpm_nonp[rowMax(as.matrix(tpm_nonp)) > 20, ]
tpm_nonp$W <- rowMeans(tpm_nonp[grep("^W", names(tpm_nonp))])
tpm_nonp$M <- rowMeans(tpm_nonp[grep("^M", names(tpm_nonp))])
tpm_nonp$log2fold <- with(tpm_nonp, log2(M + 1) - log2(W + 1))

grp <- factor(gsub("(M|W).*", "\\1", colnames(tpm_nonp)), levels = c("W", "M"))
ttests <- rowttests(as.matrix(tpm_nonp), grp)

tpm_nonp <- cbind(tpm_nonp, ttests)

tpm_nonp["Hspa2", ]
(gene_nonp <- tpm_nonp[tpm_nonp$log2fold > 0.2 & tpm_nonp$p.value < 0.2, ])

new <- list(tpm = tpm, gene_nonp = gene_nonp)
save(new, file = "temp/new.rdt")

# GO, KEGG, sex chromosome info
source("../../X/function.R")
gk <- myGK(rownames(gene_nonp))

sapply(gk$GO, function(x) x$Term[1:10])

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("chromosome_name", "external_gene_name", "description"), "external_gene_name", rownames(gene_nonp), mart)
split(biomart, biomart$chromosome_name)

# GO: Spermatogenesis
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
spermatogenesis_go <- read.delim("data/spermatogenesis.txt", header = F, stringsAsFactors = F)
spermatogenesis_genes <- getBM(attributes= "external_gene_name", "go_id", spermatogenesis_go$V1, mouse)

# GO: Transcription factors
tfs_go <- read.delim("data/tfs.txt", header = F, stringsAsFactors = F)
tfs_genes <- getBM(attributes= "external_gene_name", "go_id", tfs_go$V1, mouse)

intersect(rownames(gene_nonp), spermatogenesis_genes$external_gene_name)
intersect(rownames(gene_nonp), tfs_genes$external_gene_name)

norm1 <- within(tpm + 1, {  # Poly vs Non-poly
  RW1 = W1PLM/W1NONP; RW2 = W2PLM/W2NONP; RW3 = W3PLM/W3NONP; RM1 = M1PLM/M1NONP; RM2 = M2PLM/M2NONP; RM3 = M3PLM/M3NONP
}) %>% dplyr::select(contains("R"))

grp <- factor(gsub("R(M|W).*", "\\1", names(norm1)), levels = c("W", "M"))
ttests <- rowttests(as.matrix(norm1), grp)
ttests[ttests$p.value < 0.05, ]

# Cell types enrichment
bg <- rownames(tpm)
cells <- lapply(1:8, function(x) read.xlsx("data/TableS8RB.xlsx", stringsAsFactors = F, sheetIndex = x))

Spgonia <- read.delim("data/Spgonia", stringsAsFactors = F)
PL <- read.delim("data/PL", stringsAsFactors = F)
EL <- read.delim("data/EL", stringsAsFactors = F)
LL_Z <- read.delim("data/LL_Z", stringsAsFactors = F)
LL_Z_EP <- read.delim("data/LL_Z_EP", stringsAsFactors = F)
EP <- read.delim("data/EP", stringsAsFactors = F)
EP_LP_D <- read.delim("data/EP_LP_D", stringsAsFactors = F)
LP_D <- read.delim("data/LP_D", stringsAsFactors = F)

cells <- list(Spgonia = Spgonia, EL = EL, LL_Z = LL_Z, LL_Z_EP = LL_Z_EP, EP = EP, EP_LP_D = EP_LP_D, LP_D = LP_D)
sapply(cells, nrow)

myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bg, g2)), length(g1))
}  # Pr(count >= length(intersect(g1, g2)))

sapply(cells, function(x) myhyper(rownames(gene_nonp), x$gene_name))
lapply(cells, function(x) intersect(rownames(gene_nonp), x$gene_name))

# IN and polysomic with loose criteria and intersections
