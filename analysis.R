library(biomaRt)
library(dplyr)
library(ggvis)
library(xlsx)
library(Gviz)
library(Pviz)
library(Rsamtools)
library(seqinr)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")

load("data/ge.rdt")
ge_tpm <- sapply(ge, function(x) x$TPM) %>% as.data.frame
rownames(ge_tpm) = ge[[1]]$gene_id

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), 
  filters = "ensembl_transcript_id", values = ge[[1]]$gene_id, mart)
biomart <- biomart[match(ge[[1]]$gene_id, biomart$ensembl_transcript_id), ]

listAttributes(mart)
test <- getBM(c("peptide", "family"), filters = "external_gene_name", values = "Eif4g3", mart)

ge_tpm <- apply(ge_tpm, 2, function(x) tapply(x, biomart$external_gene_name, sum))
ge_tpm <- ge_tpm[-1, ] %>% as.data.frame

save(ge_tpm, file = "data/ge_tpm.rdt")

load("data/ge_tpm.rdt")

slp1 <- colnames(ge_tpm)
grp1 <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
spInf <- data.frame(sample = slp1, grp = factor(gsub("[123]", "", slp1), levels = grp1))
spInf <- mutate(spInf, geno = factor(gsub("^(M|W).*", "\\1", grp), levels = c("W", "M")))

ggvis1 <- function(x) {
  spInf %>% mutate(tpm = c(as.matrix(ge_tpm[x, ]))) %>%
  ggvis(~as.numeric(grp), ~tpm) %>% layer_boxplots(fill=~grp, width = 0.5) %>% layer_text(text:=~sample)
}

ggvis1("Hspa2")

tpm <- ge_tpm[apply(ge_tpm[grep("[123]IN", slp1)], 1, function(x) max(x) > 10), ]

tpm_poly <- tpm %>% select(contains("PLM"))

grp3 = factor(gsub("(M|W).*", "\\1", colnames(tpm_poly)), levels = c("W", "M"))
fit <- apply(log2(tpm_poly + 1), 1, function(x) lm(x ~ grp3))

fit.r2 <- sapply(fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(fit, function (x) summary(x)$fstatistic)
fit.et <- sapply(fit, function (x) summary(x)$coefficients[2, 1])
fit.pv <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 
fit.qv <- p.adjust(fit.pv, method = "fdr")

sig_up <- fit[fit.pv < 0.05 & fit.et > 0]
sig_dn <- fit[fit.pv < 0.05 & fit.et < 0]

names(sig_up)

source("../../X/function.R")

up_gk <- myGK(names(sig_up))
dn_gk <- myGK(names(sig_dn))

norm1 <- within(tpm + 1, {
  RW1 = W1PLM/W1NONP; RW2 = W2PLM/W2NONP; RW3 = W3PLM/W3NONP; RM1 = M1PLM/M1NONP; RM2 = M2PLM/M2NONP; RM3 = M3PLM/M3NONP
}) %>% select(contains("R"))

norm2 <- within(tpm + 1, {
  RW1 = W1PLM/W1IN; RW2 = W2PLM/W2IN; RW3 = W3PLM/W3IN; RM1 = M1PLM/M1IN; RM2 = M2PLM/M2IN; RM3 = M3PLM/M3IN
}) %>% select(contains("R"))

slp2 <- colnames(norm1)
grp2 = factor(gsub(".*(M|W).*", "\\1", slp2), levels = c("W", "M"))
spInf2 <- data.frame(sample = slp2, grp = grp2)

norm1[apply(norm1, 1, function(x) any(is.infinite(x))), ]

fit <- apply(log2(norm1 + 1), 1, function(x) lm(x ~ grp2))

fit.r2 <- sapply(fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(fit, function (x) summary(x)$fstatistic)
fit.et <- sapply(fit, function (x) summary(x)$coefficients[2, 1])
fit.pv <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 
fit.qv <- p.adjust(fit.pv, method = "fdr")

sig_up <- fit[fit.pv < 0.05 & fit.et > 0]
sig_dn <- fit[fit.pv < 0.05 & fit.et < 0]

load("../X/summary.rdt")
table %>% filter(query %in% names(sig_up))
table %>% filter(query %in% names(sig_dn))

source("../X/function.R")
up_gk <- myGK(names(sig_up))
dn_gk <- myGK(names(sig_dn))

ggvis2 <- function(x) {
  spInf2 %>% mutate(value = c(as.matrix(norm1[x, ]))) %>%
  ggvis(~as.numeric(grp), ~value) %>% layer_boxplots(fill=~grp, width = 0.5) %>% layer_text(text:=~sample)
}

ggvis2("Hspa2")

load("data/ge.rdt")
txInf <- ge[[1]] %>% select(gene_id, length)
count <- sapply(ge, function(x) x$expected_count) %>% as.data.frame
rownames(count) <- txInf$gene_id
txEif <- read.delim("shiny/gviz/Eif4g3.txt", header = F, stringsAsFactors = F)
count <- count[txEif$V1, ]
count[grep("IN", colnames(count))]
count <- count[apply(count, 1, max) > 20, ] 
(txEif <- rownames(count))  # Expressed Eif4g3 transcripts 

setwd("shiny")  # Eif4g3 DNA sequences

name <- list.files(path = "gviz/bam/", pattern = "*ALL.bam$")
lapply(name, function(x) indexBam(paste("gviz/bam", x, sep = "/")))
sample <- lapply(name, function(x) AlignmentsTrack(paste("gviz/bam", x, sep = "/")))
names(sample) <- gsub("_.*", "", name)
seqs <- SequenceTrack("gviz/seqs.fa", name = "C3H")

gvizList <- list()
gvizList$seqs <- seqs
gvizList$sample <- sample
save(gvizList, file = "gvizList.rdt")

load("gvizList.rdt")

axis <- GenomeAxisTrack()
seqs <- gvizList$seqs
sample <- gvizList$sample
plotTracks(list(seqs, axis, sample[[1]]), chromosome = seqnames(seqs)[1], from = 1, to = 2670)

dna.c3 <- read.fasta(file = "gviz/seqs.fa")
dna.mt <- read.fasta(file = "gviz/MIN.fa")

dna.diff <- sapply(txEif, function(i) {
  pos = which(dna.c3[[i]] != dna.mt[[i]])
  rbind(pos, REF = dna.c3[[i]][pos], ALT = dna.mt[[i]][pos]) %>% as.data.frame
})

pep.c3 <- getTrans(dna.c3); names(pep.c3) <- txEif
pep.mt <- getTrans(dna.mt); names(pep.mt) <- txEif

pep.diff <- sapply(txEif, function(i) {
  pos = which(pep.c3[[i]] != pep.mt[[i]])
  rbind(pos, REF = pep.c3[[i]][pos], ALT = pep.mt[[i]][pos]) %>% as.data.frame
})

pviz.axis <- ProteinAxisTrack(addNC = TRUE, littleTicks = TRUE)

pviz.seq1 <- ProteinSequenceTrack(pep.c3[[1]][1701:1750], name = "CH3")
pviz.seq2 <- ProteinSequenceTrack(pep.mt[[1]][1701:1750], name = "MT")
plotTracks(trackList = c(pviz.seq1, pviz.seq2), from = 1, to = 53)
plotTracks(trackList = c(pviz.axis, pviz.seq1, pviz.seq2), from = 1, to = 53)

setwd("../RBP")

txEif <- read.delim("../shiny/gviz/Eif4g3.txt", header = F, stringsAsFactors = F)$V1[c(2:5)]
eif_c3 <- lapply(txEif, function(x) read.xlsx("dataCSV1722_Eif4g3_C3.xlsx", stringsAsFactors = F, sheetName = x, startRow = 2))
eif_mt <- lapply(txEif, function(x) read.xlsx("dataCSV1722_Eif4g3_MT.xlsx", stringsAsFactors = F, sheetName = x, startRow = 2))
names(eif_c3) <- names(eif_mt) <- txEif
save(eif_c3, eif_mt, file = "eif.rdt")

load("eif.rdt")

eif_c3 <- lapply(eif_c3, function(x) mutate(x, UID = paste(From, To, Sequence, sep = "_")))
eif_mt <- lapply(eif_mt, function(x) mutate(x, UID = paste(From, To, Sequence, sep = "_")))

for (i in 1:4) {
  diff1 <- setdiff(eif_c3[[i]]$UID, eif_mt[[i]]$UID)
  diff2 <- setdiff(eif_mt[[i]]$UID, eif_c3[[i]]$UID); 
  eif_c3[[i]] <- filter(eif_c3[[i]], UID %in% diff1)
  eif_mt[[i]] <- filter(eif_mt[[i]], UID %in% diff2)
}

c3_2 <- lapply(eif_c3, function(x) tapply(x$Name, x$UID, function(z) paste(z, collapse = ",")) %>% as.data.frame)
mt_2 <- lapply(eif_mt, function(y) tapply(y$Name, y$UID, function(z) paste(z, collapse = ",")) %>% as.data.frame)

data(geneModels)
grtrack <- GeneRegionTrack(geneModels, chromosome = "chr7", name = "foo")
plotTracks(grtrack, transcriptAnnotation = "symbol")

x = eif_c3[[4]]
y = eif_mt[[4]]

aTrack <- AnnotationTrack(start = x$From, width = x$To - x$From, id = x$Name, name = "C3", featureAnnotation = "id", fontcolor.feature = "red")
plotTracks(aTrack, from = 300, to = 5500)
aTrack <- AnnotationTrack(start = y$From, width = y$To - y$From, id = y$Name, name = "MT", featureAnnotation = "id", fontcolor.feature = "red")
plotTracks(aTrack, from = 300, to = 5500)

genes <- c(names(sig_up), names(sig_dn))

Hspa2 <- read.delim("dataCSV8340_Hspa2.csv", stringsAsFactors = F, sep = ",", skip = 1)
Hspa2 <- Hspa2 %>% mutate(Range = paste(From, To, sep = "_"))

unique(Hspa2$Name)
table(unique(Hspa2$Name) %in% genes)
table(unique(Hspa2$Name) %in% names(fit))
unique(Hspa2$Range)
split(Hspa2$Name, Hspa2$Range)
