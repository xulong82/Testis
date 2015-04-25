library(Gviz)
library(Rsamtools)
library(mygene)

setwd("~/Dropbox/GitHub/Testis/")

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

setwd("shiny")

name <- list.files(path = "bam/", pattern = "*ALL.bam$")
lapply(name, function(x) indexBam(paste("bam", x, sep = "/")))
sample <- lapply(name, function(x) AlignmentsTrack(paste("bam", x, sep = "/")))
names(sample) <- gsub("_.*", "", name)

seqs <- SequenceTrack("../gviz/seqs.fa", name = "C3H")

seqnames(seqs)
names(sample)

alignList <- list()
alignList$seqs <- seqs
alignList$sample <- sample
alignList$tx_map <- tx_map
save(alignList, file = "alignList.rdt")

read1 <- sample[[1]]
axis <- GenomeAxisTrack()
plotTracks(list(seqs, axis, read1), chromosome = "ENSMUST00000080449", from = 1, to = 2670, add53 = T)
