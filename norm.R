# Copyright: Xulong Wang (xulong.wang@jax.org)

library(xlsx)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/FY")
load("tpm.rdt")

# --- Complete dataset with genes have TPM > 5 in at least 2 samples
data <- data.frame(row.names = tpm.dt$Symbol, tpm.dt[2:ncol(tpm.dt)])
data <- data[apply(data, 1, function (x) length(x[x>5]) > 2), ]

# --- Separate the ERCC spike ins from the others
genes <- data[! grepl("ERCC", rownames(data)), ]
spike <- data[grep("ERCC", rownames(data)), ]

save(genes, file = "Shiny/raw.rdt")
write.xlsx(genes, file = "data.xlsx", sheetName = "RAW", append = T)

# --- ERCC mol concentration vs TPM plot
mol <- read.delim("cms_095046.txt", stringsAsFactors = F)
mol <- mol[match(rownames(spike), mol$ERCC.ID), 4] * 0.02

par(mfrow = c(3, 6))
for (i in 1:ncol(spike)) plot(log2(mol), log2(spike[, i]), main = colnames(spike)[i])

# --- Normalization among samples with scaling factors from modeling the ERCC transcripts with a linear model
sf <- apply(spike[, 2:ncol(spike)], 2, function (x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 1])
sf <- c(M1IN = 1, sf)
genes <- t(apply(genes, 1, "/", sf))

save(genes, file = "Shiny/norm.rdt")
write.xlsx(genes, file = "data.xlsx", sheetName = "NEW", append = T)

# --- Relative to the IN samples
uid <- gsub("[123]", "", colnames(genes))
geno <- gsub("^(W|M).*", "\\1", uid)
cond <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 

fc <- NULL
for (idx in cond) fc <- cbind(fc, rowMeans(genes[, uid == idx]))
colnames(fc) <- cond

fc <- cbind(fc, WNONP2WIN = fc[, "WNONP"] / fc[, "WIN"], MNONP2MIN = fc[, "MNONP"] / fc[, "MIN"],
                WPLM2WIN = fc[, "WPLM"] / fc[, "WIN"], MPLM2MIN = fc[, "MPLM"] / fc[, "MIN"])

write.xlsx(as.data.frame(fc), file = "data.xlsx", sheetName = "FOLD_CHANGE", append = T)
write.table(fc, file = "temp.txt", quote = F)

# --- Gene summary
library(mygene)
geneId <- rownames(genes)
mus2hg <- read.delim("~/Dropbox/X/hg2mus.map", header = F, stringsAsFactors = F)
table <- queryMany(geneId, scopes="symbol", species="mouse", fields = c("name", "summary"))
table <- table[!duplicated(table$query), ]
table$hg <- mus2hg$V1[match(geneId, mus2hg$V4)]
query.hg <- queryMany(table$hg, scopes="symbol", species="human", fields = c("name", "summary"))
query.hg <- query.hg[!duplicated(query.hg$query), ]
summary_hg <- query.hg$summary[match(table$hg, query.hg$query)]
table$summary <- paste("Mus:", table$summary, "Hg:", summary_hg)
table <- table[c("query", "name", "hg", "summary")]
rownames(table) <- NULL
table <- as.data.frame(table)
save(table, file = "./Shiny/summary.rdt")

# # --- MAP to mol concentration 
# genes <- log2(genes + 1)
# for (i in 1:ncol(genes)) {
#   coef <- summary(lm(log2(spike[, i]) ~ log2(mol)))$coefficients
#   genes[, i] <- (genes[, i] - coef[1, 1]) / coef[2, 1]
# }
# gene1 <- genes
# genes <- 2^genes
# 
# # --- SCALE in samples ---
# gene <- apply(gene, 2, scale)
# rownames(gene1) <- rownames(genes)
# genes <- gene1
# 
# # --- RUVSeq method ---
# library(RUVSeq)
# library(EDASeq)
# x <- as.factor(gsub("[123]", "", colnames(genes)))
# set <- newSeqExpressionSet(as.matrix(data), phenoData = data.frame(x, row.names=colnames(data)))
# plotRLE(set)
# set <- betweenLaneNormalization(set, which = "upper")
# plotRLE(set)
# set <- RUVg(set, rownames(spike), k = 1)
# plotRLE(set)
# data <- normCounts(set)
