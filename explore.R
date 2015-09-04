library(dplyr)
library(biomaRt)
library(preprocessCore)
library(xlsx)
library(genefilter)
library(Biobase)
library(matrixStats)
library(limma)
library(igraph)

# Q: what leads Eif4g3 mutation to infertility?
rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")
source("../../X/function.R")
load("data/dataList.rdt")
for(obj in names(dataList)) assign(obj, dataList[[obj]])
load("data/ge.rdt")
tpm <- sapply(ge, function(x) x$TPM) %>% as.data.frame
tx_select <- ge[[1]]$gene_id[rowMaxs(as.matrix(tpm)) > 20] # pool in transcript

data <- tpm[! grepl("ERCC", rownames(tpm)), ]
gene <- data[biomart$ensembl_transcript_id, ]
gene_exclude <- data[! rownames(data) %in% biomart$ensembl_transcript_id, ]
gene <- apply(gene, 2, function(x) tapply(x, biomart$external_gene_name, sum))
bg <- rownames(gene)[rowMax(gene) > 30] # pool in gene

grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
grp2 <- factor(gsub("[123]", "", colnames(gene)), levels = grp1)
matboxplot(gene, groupFactor = grp2)
gene <- gene[! rowMaxs(gene) > 1e4, ]  # removal: tpm
gene <- sweep(gene, 2, colSums(gene), "/") * 1e6 # tpm
boxplot(as.matrix(gene)["Hspa2", ] ~ grp2, col = 2:7)

select <- "IN"
select <- "PLM"
select <- "NONP"

gene_select <- as.data.frame(gene) %>% dplyr::select(contains(select))  
gene_select <- gene_select[rowMax(as.matrix(gene_select)) > 20, ] # tpm
gene_select$W <- rowMeans(gene_select[grep("^W", names(gene_select))])
gene_select$M <- rowMeans(gene_select[grep("^M", names(gene_select))])
gene_select$log2fold <- with(gene_select, log2(M + 1) - log2(W + 1))
group <- factor(gsub("(M|W).*", "\\1", names(gene_select[1:6])), levels = c("W", "M"))
(design <- model.matrix(~ group))
v <- voom(gene_select[1:6], design, plot = T)
ebayesFit <- lmFit(v, design) %>% eBayes
gene_select$p.value <- ebayesFit$p.value[, "groupM"]
gene_select$q.value <- qvalue(gene_select$p.value)$qvalues

gene_select <- gene_select[gene_select$p.value < 0.05, ]
gene_select <- gene_select[gene_select$p.value < 0.05 & gene_select$log2fold > 0.2, ]
gene_select <- gene_select[gene_select$p.value < 0.05 & gene_select$log2fold < -0.2, ]
gene_select <- gene_select[order(gene_select$q.value), ]

input = rownames(gene_select)
gene_select_gk <- myGK(input)
sapply(gene_select_gk$GO, function(x) x$Term[1:20]) # top terms
sapply(gene_select_gk$GO[c(1, 3)], function(x) x$Term[1:20]) # top terms

boxplot(as.matrix(gene_select)["Actr10", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["Dynlt1b", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["Cep70", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["Gm3376", 1:6] ~ group, col = 3:2)

write.xlsx(gene_select, file = "new.xlsx", sheetName = "gene_select", append = T)
write.xlsx(gene_select_gk$GO$BP, file = "new.xlsx", sheetName = "gene_select_GO_BP", append = T)
write.xlsx(gene_select_gk$GO$MF, file = "new.xlsx", sheetName = "gene_select_GO_MF", append = T)
write.xlsx(gene_select_gk$GO$CC, file = "new.xlsx", sheetName = "gene_select_GO_CC", append = T)
write.xlsx(gene_select_gk$KEGG, file = "new.xlsx", sheetName = "gene_select_KEGG", append = T)

file <- "pathway/IN_new.tsv" # Network iRegulon::IN genes 
ireg <- read.delim(file, comment.char = ";", stringsAsFactors = F)
factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))
factor <- lapply(factor, function(x) x[x %in% bg])
target <- target[sapply(factor, length) > 0]
factor <- factor[sapply(factor, length) > 0]
edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
edges <- do.call(rbind, edges)
edges <- edges[! duplicated(edges), ]
igraph.dt <- graph.data.frame(edges)  # IGRAPH
igraph.dt$layout <- layout.sphere
V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% unlist(factor)] <- "gold"
V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.5
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"

pdf(file = "IN_new.pdf", width = 12, height = 12)
plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)
dev.off()

# Eif4g3 & Hspa2 protein interaction
pbp_Eif4g3 <- read.delim("String/tabdelimited.Eif4g3.txt", stringsAsFactors = F)
pbp_Hspa2 <- read.delim("String/tabdelimited.Hspa2.txt", stringsAsFactors = F)

intersect(input, unique(c(as.matrix(pbp_Eif4g3[, 1:2]))))
x = intersect(input, unique(c(as.matrix(pbp_Hspa2[, 1:2])))) # 6

pdf("Pdf/Hspa2_pbp.pdf", width = 8, height = 5)
par(mfrow = c(2, 3))
lapply(x, function(x) boxplot(as.matrix(gene_select)[x, 1:6] ~ group, col = 3:2, main = x))
dev.off()

# Pattern of NONP genes in development
load("meiotic/Y_log2TPM.RData")
attributes <- c("ensembl_gene_id", "external_gene_name")
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(attributes, "ensembl_gene_id", rownames(Y), ensembl)
meiotic <- Y[biomart$ensembl_gene_id, ]
meiotic <- apply(meiotic, 2, function(x) tapply(x, biomart$external_gene_name, sum))
meiotic <- meiotic[rowMax(meiotic) > 5, ] # 8063 genes

mei_group <- gsub("_.*", "", colnames(meiotic)) 
mei_stage <- unique(mei_group)
mei_group_fac <- factor(mei_group, levels = mei_stage)
meiotic_mean <- sapply(mei_stage, function(x) rowMeans(meiotic[, mei_group == x]))

pca.ge <- meiotic
pca.ge <- pca.ge - rowMeans(pca.ge[, grep("^8_", colnames(pca.ge))])
pca <- prcomp(pca.ge)

input_select <- intersect(rownames(meiotic), input)
pca.ge_select <- meiotic[input_select, ]
pca.ge_select <- pca.ge_select - rowMeans(pca.ge_select[, grep("^8_", colnames(pca.ge_select))])
pca_select <- prcomp(pca.ge_select)

barplot(pca$sdev / sum(pca$sdev), col = c(rep("red", 6), rep("grey", 29)), xlab = "PC", ylab = "% variants explained")
barplot(pca_select$sdev / sum(pca_select$sdev), col = c(rep("red", 6), rep("grey", 29)), xlab = "PC", ylab = "% variants explained")

pdf("Pdf/pca1.pdf", width = 8, height = 5)
par(mfrow = c(2, 3), mar = c(2, 4, 4, 2))
lapply(1:6, function(x) barplot(pca_select$rotation[, x], col = as.numeric(mei_group_fac), xaxt = "n", main = x))
dev.off()
  
# strongest signal is genes that keep increasing!!!
# this is puzzling!!! how to understand? Find a gene group that don't show same dominant signal

pca_gene_pc <- pca_select$x[, 1]
pca_gene_pc <- pca_select$x[, 2]

yy1 <- which(pca_gene_pc > sd(pca_gene_pc)) %>% names
yy2 <- which(pca_gene_pc < -sd(pca_gene_pc)) %>% names

graph.dt <- meiotic_mean[yy1, ]
graph.dt <- meiotic_mean[yy2, ]
graph.dt <- meiotic_mean[input_select, ]
graph.dt <- graph.dt - graph.dt[, 1]

pdf("Pdf/pca3.pdf", width = 8, height = 5)
plot(1:6, rep(0, 6), type = "l", col = "blue", xlab = "", xaxt = "n", ylab = "", ylim = c(-2, 9))
plot(1:6, rep(0, 6), type = "l", col = "blue", xlab = "", xaxt = "n", ylab = "", ylim = c(-2, 2))
lapply(1:nrow(graph.dt), function(i) lines(graph.dt[i, ], type = "l", col = "grey30"))
lines(colMeans(graph.dt), type = "b", lwd = 3, col = "red")
axis(1, at=1:6, labels = paste0(c(8, 10, 12, 14, 16, 18), "d"))
dev.off()

# shared 3'UTR RBP
genes <- bg # to exact RNAmap database
genes <- rownames(gene_select)
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attributes <- c("external_gene_name", "ensembl_transcript_id", "chromosome_name", "strand", "3_utr_start", "3_utr_end")
mouse_mart <- getBM(attributes, "external_gene_name", genes, mouse)
mouse_mart <- filter(mouse_mart, ensembl_transcript_id %in% tx_select)
mouse_mart <- filter(mouse_mart, ! (is.na(mouse_mart$`3_utr_start`) | is.na(mouse_mart$`3_utr_end`)))
mouse_mart$strand_2 <- "+"
mouse_mart$strand_2[mouse_mart$strand == -1] = "-"

locationAll <- with(mouse_mart, paste0("chr", chromosome_name, ":", `3_utr_start`, "-", `3_utr_end`, ":", strand_2))
write(locationAll[1:4956], file = "RBPmap/locationAll_part1.txt") # RBPmap
write(locationAll[4957:9913], file = "RBPmap/locationAll_part2.txt") # RBPmap

location_select <- with(mouse_mart, paste0("chr", chromosome_name, ":", `3_utr_start`, "-", `3_utr_end`, ":", strand_2))
write(location_select, file = "RBPmap/location.txt") # RBPmap

rbp_Hspa2 <- read.delim("RBPmap/RBP_Hspa2.txt", header = F, stringsAsFactors = F)$V1
rbp_All <- read.delim("RBPmap/RBP_All.txt", header = F, stringsAsFactors = F)$V1

rbp_Hspa2 <- read.delim("RBPmap/RBP_Hspa2_stringent.txt", header = F, stringsAsFactors = F)$V1
rbp_All <- read.delim("RBPmap/RBP_All_stringent.txt", header = F, stringsAsFactors = F)$V1

sort(sapply(rbp_Hspa2, function(x) sum(rbp_All == x)), decreasing = T)

tx_select_3p <- mouse_mart$ensembl_transcript_id
utr3p <- getSequence(seqType='3utr', mart = mouse, type="ensembl_transcript_id", id = tx_select_3p)
utr3p <- utr3p[utr3p$`3utr` != "Sequence unavailable", ]
utr3p <- utr3p[! grepl("No UTR", utr3p$`3utr`), ]
write(c(t(utr3p[, c(2, 1)])), file = "RBP/utr3p_nonpoly_new.txt") # CISBP

# RAID: RNA-associated interaction database. Symbol are confusing!!! Author are correcting.
list.files("RAID")  
rna_rna <- read.xlsx("RAID/RNA-RNA-download.xlsx", sheetIndex = 1, stringsAsFactors = F)
rna_protein <- read.xlsx("RAID/RNA-Protein-download.xlsx", sheetIndex = 1, stringsAsFactors = F)
dataList$RAID$rna_rna <- rna_rna
dataList$RAID$rna_protein <- rna_protein
