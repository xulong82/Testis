library(biomaRt)
library(ggplot2)
library(igraph)
library(dplyr)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")
load("../X/summary.rdt")
summary <- function(x) table %>% filter(query %in% x)

# LOAD DATA
load("data/newList.rdt")
# for(obj in names(newList)) assign(obj, newList[[obj]])
for(obj in names(newList$more)) assign(obj, newList$more[[obj]])
load("data/biology.rdt")
for(obj in names(biology)) assign(obj, biology[[obj]])
geneId1 <- rownames(tpm)[apply(tpm, 1, function(x) max(x) > 50)]

load("data/ge.rdt")
ge_tpm <- sapply(ge, function(x) x$TPM) %>% as.data.frame
rownames(ge_tpm) = ge[[1]]$gene_id
txId1 <- rownames(ge_tpm)[apply(ge_tpm, 1, function(x) max(x) > 20)]

lapply(nonpoly, function(x) names(x))
lapply(nonpoly, function(x) summary(names(x)))

# Eif4g3 RBP in the RNA level
rbp$Eif4g3$c3 <- lapply(rbp$Eif4g3$c3, function(x) mutate(x, UID = paste(From, To, Sequence, sep = "_")))
rbp$Eif4g3$mt <- lapply(rbp$Eif4g3$mt, function(x) mutate(x, UID = paste(From, To, Sequence, sep = "_")))

rbp$Eif4g3$c3 <- lapply(rbp$Eif4g3$c3, function(x) filter(x, Score > 10))
rbp$Eif4g3$mt <- lapply(rbp$Eif4g3$mt, function(x) filter(x, Score > 10))

for (i in 1:4) {
  diff1 <- setdiff(rbp$Eif4g3$c3[[i]]$UID, rbp$Eif4g3$mt[[i]]$UID)
  diff2 <- setdiff(rbp$Eif4g3$mt[[i]]$UID, rbp$Eif4g3$c3[[i]]$UID); 
  rbp$Eif4g3$diff$c3[[i]] <- filter(rbp$Eif4g3$c3[[i]], UID %in% diff1)
  rbp$Eif4g3$diff$mt[[i]] <- filter(rbp$Eif4g3$mt[[i]], UID %in% diff2)
}

names(rbp$Eif4g3$diff$c3) <- names(rbp$Eif4g3$diff$mt) <- names(rbp$Eif4g3$c3)

c3Rbp <- do.call(rbind, rbp$Eif4g3$diff$c3)$Name %>% unique
mtRbp <- do.call(rbind, rbp$Eif4g3$diff$mt)$Name %>% unique
intersect(c3Rbp, mtRbp)
Eif4g3Rbp <- union(c3Rbp, mtRbp)
(Eif4g3Rbp <- intersect(Eif4g3Rbp, geneId1))

lapply(poly, function(x) intersect(Eif4g3Rbp, names(x)))
lapply(nonpoly, function(x) intersect(Eif4g3Rbp, names(x)))

# Hspa2 RBP in the RNA level
# rbp$Hspa2 <- read.delim("RBP/dataCSV8340_Hspa2.csv", stringsAsFactors = F, sep = ",", skip = 1)
rbp$Hspa2 <- rbp$Hspa2 %>% mutate(Range = paste(From, To, sep = "_"))
rbp$Hspa2 <- filter(rbp$Hspa2, Score > 10)
Hspa2Rbp <- unique(rbp$Hspa2$Name)
(Hspa2Rbp <- intersect(Hspa2Rbp, geneId1))
split(rbp$Hspa2$Name, rbp$Hspa2$Range)

intersect(Eif4g3Rbp, Hspa2Rbp)  # No specificity. Potentially interesting! Not likely!

lapply(poly, function(x) intersect(Hspa2Rbp, names(x)))  # Potentially interesting!
lapply(nonpoly, function(x) intersect(Hspa2Rbp, names(x)))
lapply(In, function(x) intersect(Hspa2Rbp, names(x)))

# RBP in 3'UTR
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
sy2tx <- getBM(c("external_gene_name", "ensembl_transcript_id"), "external_gene_name", names(nonpoly$up), mouse)
sy2tx <- getBM(c("external_gene_name", "ensembl_transcript_id"), "external_gene_name", names(norm$down), mouse)

utr3p <- getSequence(seqType='3utr', mart = mouse, type="ensembl_transcript_id", id = sy2tx[, 2])
utr3p <- utr3p[utr3p[, 1] != "Sequence unavailable", ]
utr3p <- utr3p[utr3p[, 2] %in% txId1, ]

write(c(t(utr3p[, c(2, 1)])), file = "RBP/utr3p_nonpoly.txt")
write(c(t(utr3p[, c(2, 1)])), file = "RBP/utr3p_norm.txt")

utr3p_rbp <- lapply(1:11, function(x) read.xlsx("RBP/dataCSV104_Nonpoly_UTR3.xlsx", sheetIndex = x, startRow = 2, stringsAsFactors = F))
names(utr3p_rbp) <- utr3p[, 2]
newList$utr3p$nonpoly <- utr3p_rbp

utr3p_rbp <- lapply(1:9, function(x) read.xlsx("RBP/dataCSV104_Norm_UTR3.xlsx", sheetIndex = x, startRow = 2, stringsAsFactors = F))
names(utr3p_rbp) <- utr3p[, 2]
newList$utr3p$norm <- utr3p_rbp

utr3p_rbp <- lapply(utr3p$nonpoly, function(x) filter(x, Score > 10))
utr3p_rbp <- lapply(utr3p$norm, function(x) filter(x, Score > 10))

utr3p_rbp_hspa2 <- utr3p_rbp[["ENSMUST00000080449"]]
utr3p_rbp_hspa2 <- unique(utr3p_rbp_hspa2$Name)
(utr3p_rbp_hspa2 <- intersect(utr3p_rbp_hspa2, geneId1))

lapply(poly, function(x) intersect(utr3p_rbp_hspa2, names(x)))  # Potentially interesting!
lapply(nonpoly, function(x) intersect(utr3p_rbp_hspa2, names(x)))
lapply(In, function(x) intersect(utr3p_rbp_hspa2, names(x)))

lapply(utr3p_rbp_hspa2, function(x) sapply(utr3p_rbp, function(y) x %in% y$Name))
sort(sapply(utr3p_rbp_hspa2, function(x) sum(sapply(utr3p_rbp, function(y) x %in% y$Name))), decreasing = T)

# Eif4g3 & Hspa2 protein interaction
# biology$pbp$Eif4g3 <- read.delim("String/tabdelimited.Eif4g3.txt", stringsAsFactors = F)
# biology$pbp$Hspa2 <- read.delim("String/tabdelimited.Hspa2.txt", stringsAsFactors = F)

lapply(poly, function(x) intersect(names(x), unique(c(as.matrix(pbp$Eif4g3[, 1:2])))))  # Potentially interesting
lapply(nonpoly, function(x) intersect(names(x), unique(c(as.matrix(pbp$Eif4g3[, 1:2])))))
lapply(norm, function(x) intersect(names(x), unique(c(as.matrix(pbp$Eif4g3[, 1:2])))))

lapply(poly, function(x) intersect(names(x), unique(c(as.matrix(pbp$Hspa2[, 1:2])))))  # Potentially interesting
lapply(nonpoly, function(x) intersect(names(x), unique(c(as.matrix(pbp$Hspa2[, 1:2])))))
lapply(norm, function(x) intersect(names(x), unique(c(as.matrix(pbp$Hspa2[, 1:2])))))

intersect(unique(c(as.matrix(pbp$Eif4g3[, 1:2]))), unique(c(as.matrix(pbp$Hspa2[, 1:2]))))

intersect(unique(c(as.matrix(pbp$Eif4g3[, 1:2]))), c3Rbp)
intersect(unique(c(as.matrix(pbp$Eif4g3[, 1:2]))), mtRbp)
intersect(unique(c(as.matrix(pbp$Eif4g3[, 1:2]))), Hspa2Rbp)

lapply(rbp$Eif4g3$diff$c3, function(x) x %>% filter(Name == "Pabpc1"))
lapply(rbp$Eif4g3$diff$mt, function(x) x %>% filter(Name == "Pabpc1"))

# GO: Spermatogenesis
biology$spermatogenesis$go <- read.delim("data/spermatogenesis.txt", header = F, stringsAsFactors = F)
biology$spermatogenesis$genes <- getBM(attributes= "external_gene_name", "go_id", Eif4g3$spermatogenesis$go$V1, mouse)

# GO: Transcription factors
biology$tfs$go <- read.delim("data/tfs.txt", header = F, stringsAsFactors = F)
biology$tfs$genes <- getBM(attributes= "external_gene_name", "go_id", Eif4g3$tfs$go$V1, mouse)

goSearch = getBM(c("go_id", "name_1006", "definition_1006"), "external_gene_name", c("Hspa2", "Hspa8"), mouse)

(interPoly <- lapply(poly, function(x) intersect(spermatogenesis$genes[, 1], names(x))))
(interNonpoly <- lapply(nonpoly, function(x) intersect(spermatogenesis$genes[, 1], names(x))))
lapply(norm, function(x) intersect(spermatogenesis$genes[, 1], names(x)))
summary(interPoly$up)
summary(interPoly$down)

# Eif4g & Hspa gene families
(Hspa <- rownames(tpm)[grep("Hspa", rownames(tpm))])
(Eif4g <- rownames(tpm)[grep("Eif4g", rownames(tpm))])

grp <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
spInf <- data.frame(sample = colnames(tpm), grp = factor(gsub("[123]", "", colnames(tpm)), levels = grp))
spInf <- mutate(spInf, geno = factor(gsub("^(M|W).*", "\\1", grp), levels = c("W", "M")))

graph1 <- function(x) spInf %>% mutate(value = c(as.matrix(tpm[x, ]))) %>%
  ggplot(aes(x = grp, y = value)) + geom_boxplot(aes(fill = geno)) +
  scale_fill_manual(values = c("dodgerblue3", "firebrick1")) +
  theme_bw() + xlab("") + ylab("") + ggtitle(x)

lapply(Hspa, graph1)
lapply(Eif4g, graph1)
lapply(names(nonpoly$up), graph1)
lapply(names(norm$down), graph1)

# MASTER REGULATORS OF POLYSOMIC GENES
write(names(less$poly$up), file = "pathway/polyUp.txt")
write(names(less$poly$down), file = "pathway/polyDown.txt")
file <- "pathway/polyUp_iregulon.tsv"
file <- "pathway/polyDown_iregulon.tsv"
universe <- rownames(tpm)[apply(tpm, 1, function(x) max(x) > 20)]
ireg <- read.delim(file, comment.char = ";", stringsAsFactors = F)
factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))
factor <- lapply(factor, function(x) x[x %in% universe])
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

pdf(file = "polyUp.pdf", width = 12, height = 12)
pdf(file = "polyDown.pdf", width = 12, height = 12)
plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)
dev.off()
