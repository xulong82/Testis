library(dplyr)
library(biomaRt)
library(quantro)
library(preprocessCore)
library(xlsx)
library(genefilter)
library(Biobase)
library(matrixStats)
library(limma)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")
source("../../X/function.R")
load("data/dataList.rdt")
for(obj in names(dataList)) assign(obj, dataList[[obj]])

# GO::Spermatogenesis
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
spermatogenesis_go <- read.delim("data/spermatogenesis.txt", header = F, stringsAsFactors = F)
spermatogenesis_genes <- getBM(attributes= "external_gene_name", "go_id", spermatogenesis_go$V1, mart)

# ---
load("data/ge.rdt")
cnt <- sapply(ge, function(x) x$expected_count) %>% as.data.frame
tpm <- sapply(ge, function(x) x$TPM) %>% as.data.frame
rownames(cnt) <- rownames(tpm) <- ge[[1]]$gene_id

query <- rownames(cnt)[! grepl("ERCC", rownames(cnt))]
attributes <- c("ensembl_transcript_id", "external_gene_name")
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(attributes, "ensembl_transcript_id", query, ensembl)
# ---

data <- tpm[! grepl("ERCC", rownames(tpm)), ]
data <- cnt[! grepl("ERCC", rownames(cnt)), ]

gene <- data[biomart$ensembl_transcript_id, ]
gene_exclude <- data[! rownames(data) %in% biomart$ensembl_transcript_id, ]
gene <- apply(gene, 2, function(x) tapply(x, biomart$external_gene_name, sum))

grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
grp2 <- factor(gsub("[123]", "", colnames(gene)), levels = grp1)
matboxplot(gene, groupFactor = grp2)
colMaxs(gene)
colMaxs(gene) / colSums(gene)

gene[rowMaxs(gene) > 1e4, ]
gene[rowMaxs(gene) > 2e5, ]
gene[rowMaxs(gene) > 1e6, ]

gene <- gene[! rowMaxs(gene) > 1e4, ]  # removal: tpm
gene <- gene[! rowMaxs(gene) > 2e5, ]  # removal: count

# re-normalize
gene <- sweep(gene, 2, colSums(gene), "/") * 1e6 # tpm
gene <- sweep(gene, 2, colSums(gene), "/") * 5e7 # count
(qtest <- quantro(object = log2(tpm + 1), groupFactor = grp2, B = 1e3))  # quantro

# Non-polysomic
gene_nonp <- as.data.frame(gene) %>% dplyr::select(contains("NONP"))  # Non-poly
gene_nonp <- gene_nonp[rowMax(as.matrix(gene_nonp)) > 20, ] # tpm
gene_nonp <- gene_nonp[rowMax(as.matrix(gene_nonp)) > 300, ] # count
group <- factor(gsub("(M|W).*", "\\1", colnames(gene_nonp)), levels = c("W", "M"))

gene_nonp$W <- rowMeans(gene_nonp[grep("^W", names(gene_nonp))])
gene_nonp$M <- rowMeans(gene_nonp[grep("^M", names(gene_nonp))])
gene_nonp$log2fold <- with(gene_nonp, log2(M + 1) - log2(W + 1))

matboxplot(log2(gene_nonp[1:6] + 0.5), groupFactor = group)
qqplot(log2(gene_nonp$W + 1), log2(gene_nonp$M + 1), xlab = "WT", ylab = "MT")
abline(0, 1)

# ttests::tpm
ttests <- rowttests(as.matrix(log2(gene_nonp[1:6] + 1)), group)
gene_nonp_ttest <- cbind(gene_nonp, ttests)

gene_nonp_ttest["Hspa2", ]
ttest_select <- gene_nonp_ttest[with(gene_nonp_ttest, log2fold > 0.2 & p.value < 0.05), ]
(ttest_id <- rownames(ttest_select))
gk_ttest <- myGK(ttest_id)
sapply(gk_ttest$GO, function(x) x$Term[1:10])

gene_nonp_scale <- apply(gene_nonp, 1, scale) %>% t
gene_nonp_quantile <- normalize.quantiles(as.matrix(gene_nonp))
dimnames(gene_nonp_quantile) <- dimnames(gene_nonp_scale) <- dimnames(gene_nonp)

# ebayes::tpm/count
(design <- model.matrix(~ group))
v <- voom(gene_nonp[1:6], design, plot = T)
v <- voom(gene_nonp[1:6], design, plot = T, normalize.method = "quantile") # no reason to do so
matboxplot(v$E, groupFactor = group)

lmFit <- lmFit(v, design)
ebayesFit <- eBayes(lmFit)
(top <- topTable(ebayesFit,coef=ncol(design), n = 20))

qvalue <- qvalue(ebayesFit$p.value[, "groupM"])$qvalues
which(qvalue < 0.05)

gene_ebayes <- which(qvalue < 0.05 & lmFit$coefficients[, "groupM"] > 0.2)
gene_ebayes <- which(ebayesFit$p.value[, "groupM"] < 0.05 & lmFit$coefficients[, "groupM"] > 0.2)
ebayes_id <- names(gene_ebayes)
ebayes_id_count <- names(gene_ebayes)
ebayes_id_qvalue <- names(gene_ebayes)
ebayes_id_qnorm <- names(gene_ebayes)
lapply(ebayes_id_qvalue, function(y) boxplot(as.matrix(gene_nonp)[y, 1:6] ~ group, col = 3:2, main = y))

# use tpm or count?
(x1 <- setdiff(ebayes_id, ebayes_id_count))
gk_x1 <- myGK(x1)
(x2 <- setdiff(ebayes_id_count, ebayes_id))
gk_x2 <- myGK(x2)

sapply(gk_x1$GO, function(x) x$Term[1:10]) # many meiotic genes, argue for tpm
sapply(gk_x2$GO, function(x) x$Term[1:10]) # irrelevant terms, argue against count

# quantile.normalize necessary in voom()??? 
x1 <- setdiff(ebayes_id, ebayes_id_qnorm)
gk_x1 <- myGK(x1)
x2 <- setdiff(ebayes_id_qnorm, ebayes_id)
gk_x2 <- myGK(x2)

sapply(gk_x1$GO, function(x) x$Term[1:10]) # many meiotic genes
sapply(gk_x2$GO, function(x) x$Term[1:10]) # only 1,low level, quantile norm isn't necessary

# use ttest::tpm or ebayes::tpm
t_ebayes <- ebayesFit$t
t_ttests <- rowttests(v$E, group)
plot(t_ebayes[, "groupM"], -t_ttests[, "statistic"], xlim = c(-10, 10), ylim = c(-10, 10))
abline(0, 1)

(x1 <- setdiff(ebayes_id, ttest_id))
gk_x1 <- myGK(x1)
(x2 <- setdiff(ttest_id, ebayes_id))
gk_x2 <- myGK(x2)

sapply(gk_x1$GO, function(x) x$Term[1:10]) # many meiotic genes, argue for ebayes
sapply(gk_x2$GO, function(x) x$Term[1:10]) # irrelevant terms, argue against ttest

### Conclusion: use TPM and ebayes method
intersect(ttest_id, spermatogenesis_genes$external_gene_name) 
intersect(ebayes_id, spermatogenesis_genes$external_gene_name)
intersect(ebayes_id_qnorm, spermatogenesis_genes$external_gene_name)
intersect(ebayes_id_count, spermatogenesis_genes$external_gene_name)

# What are the targeting chromosomes?
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mart_chr <- getBM(c("chromosome_name", "external_gene_name", "description"), "external_gene_name", ebayes_id, mart)
(chrList <- split(mart_chr, mart_chr$chromosome_name))
sapply(chrList, nrow)

# why Gm10282, Supv3l1 not in ebayes list? A: low level makes moderated deviation large!!!
boxplot(as.matrix(gene_nonp)["Supv3l1", 1:6] ~ group, col = 3:2, ylim = c(0, 100))
boxplot(as.matrix(gene_nonp)["Gm10282", 1:6] ~ group, col = 3:2, ylim = c(0, 100))

boxplot(as.matrix(gene_select)["Stip1", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["Tuba3a", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["Hspa2", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_select)["1600020E01Rik", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_nonp)["Actb", 1:6] ~ group, col = 3:2)
boxplot(as.matrix(gene_nonp)["Rad1", 1:6] ~ group, col = 3:2)
