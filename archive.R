setwd("~/Dropbox/GitHub/Testis/")

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
