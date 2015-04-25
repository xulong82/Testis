rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")

packages <- c("tidyr","dplyr", "ggvis")
sapply(packages, require, character.only = T)

trim <- read.delim("data/trim.txt", header = F, stringsAsFactors = F)
rsem <- read.delim("data/rsem_c3h.txt", header = F, stringsAsFactors = F) 

qc <- full_join(trim, rsem, by = "V1")[, -3] 
colnames(qc) <- c("sample", "trim", "count", "bowtie")

group <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
qc$sample <- gsub("_.*", "", qc$sample)
qc$group <- gsub("[123]", "", qc$sample) %>% factor(levels = group)
qc$genotype <- gsub("^(M|W).*", "\\1", qc$sample) %>% factor(levels = c("W", "M"))

qc$count <- qc$count * 1e-6
qc$trim <- as.numeric(gsub("%", "", qc$trim)) * 1e-2
qc$bowtie<- as.numeric(gsub("%", "", qc$bowtie)) * 1e-2

qc$aligned <- qc$count * qc$bowtie
qc$group1 <- as.numeric(qc$group)

ggvis_boxplots <- function(x) {x %>% add_axis("x", values = NULL, title = "Group") %>%
  layer_boxplots(fill=~group, width = 0.3) %>% layer_text(text:=~sample)}

qc %>% ggvis(~group1, ~trim) %>% scale_numeric("y", domain=c(0.5, 1)) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~bowtie) %>% scale_numeric("y", domain=c(0.5, 1)) %>% ggvis_boxplots() 
qc %>% ggvis(~group1, ~count) %>% scale_numeric("y", domain=c(25, 45)) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~aligned) %>% scale_numeric("y", domain=c(15, 35)) %>% ggvis_boxplots()

qc %>% ggvis(~count, ~bowtie) %>% layer_points(size:=1e2, fill=~group) %>% layer_text(text:=~sample)

# setwd("/data/xwang/Testis")  # CADILLAC
# name <- list.files(path = "./RSEM/", pattern = "*.results")
# ge <- lapply(name, function(x) {
#   cat(x, "\n"); filepath <- file.path("./RSEM/", x)
#   read.delim(filepath, stringsAsFactors = F)
# }); names(ge) <- gsub("_.*", "", name)
# save(ge, file = "~/Dropbox/GitHub/Testis/data/ge.rdt")

load("data/ge.rdt"); load("../X/ensembl_mus.rdt")
txInf <- ge[[1]] %>% select(gene_id, length) %>% mutate(symbol = ens.map[gene_id, 2])
ge_TPM <- sapply(ge, function(x) x$TPM) %>% as.data.frame
ge_count <- sapply(ge, function(x) x$expected_count) %>% as.data.frame
rownames(ge_count) <- rownames(ge_TPM) <- txInf$gene_id

data <- ge_count; data <- data[, qc$sample] 

genes <- data[! grepl("ERCC", rownames(data)), ]
spike <- data[grep("ERCC", rownames(data)), ]

qc$f_count <- qc$count / qc$count[1]
qc$f_aligned <- qc$aligned / qc$aligned[1]

qc$f_spike1 <- colSums(spike) / colSums(spike)[1]
qc$f_spike2 <- c(W1NONP = 1, apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 1]))

txInf1 <- txInf %>% filter(! grepl("ERCC", gene_id))
raw <- apply(genes, 2, function(x) tapply(x, txInf1$symbol, sum))
raw <- raw[apply(raw, 1, function (x) length(x[x > 20]) > 2), ]

total <- t(apply(raw, 1, "/", qc$f_count))
aligned <- t(apply(raw, 1, "/", qc$f_aligned))
spike2 <- t(apply(raw, 1, "/", qc$f_spike2))

geList <- list()
geList$raw = raw; geList$total = total; geList$aligned = aligned; geList$spike2 = spike2; 
save(geList, file = "Shiny/geList.rdt")

qc$f_spike3 <- c(W1NONP = 1, apply(spike2[, -1], 2, function(x) summary(lm(x ~ -1 + spike2[, 1]))$coefficients[1, 1]))

std_err <- apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 2])
pval <- apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 4])

qc %>% ggvis(~group1, ~f_spike1) %>% scale_numeric("y", domain=c(.5, 2.0)) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~f_spike2) %>% scale_numeric("y", domain=c(.5, 2.0)) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~f_spike3) %>% scale_numeric("y", domain=c(.5, 2.0)) %>% ggvis_boxplots()

mol <- read.delim("data/cms_095046.txt", stringsAsFactors = F)
mol <- mol[match(rownames(spike), mol$ERCC.ID), 4] * 0.02

par(mfrow = c(3, 6))
for (i in 1:ncol(spike)) plot(log2(mol), log2(spike[, i]), ylim = c(2, 14), main = names(spike)[i])

ensId <- txInf %>% filter(symbol == "Actb") %>% select(gene_id)
qc$actb <- colSums(genes[ensId$gene_id, ])

qc <- mutate(qc, actb_f_count = actb / f_count, actb_f_aligned = actb / f_aligned,
             actb_f_spike1 = actb / f_spike1, actb_f_spike2 = actb / f_spike2)

qc %>% ggvis(~group1, ~actb) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~actb_f_count) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~actb_f_aligned) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~actb_f_spike1) %>% ggvis_boxplots()

ensId <- txInf %>% filter(symbol == "Hspa2") %>% select(gene_id)
qc$hspa2 = c(t(genes[ensId$gene_id, ]))

qc <- mutate(qc, hspa2_f_count = hspa2 / f_count, hspa2_f_aligned = hspa2 / f_aligned,
             hspa2_f_spike1 = hspa2 / f_spike1, hspa2_f_spike2 = hspa2 / f_spike2)

qc %>% ggvis(~group1, ~hspa2) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~hspa2_f_count) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~hspa2_f_aligned) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~hspa2_f_spike1) %>% ggvis_boxplots()

ensId <- txInf %>% filter(symbol == "Actb") %>% select(gene_id)
qc$actb2 <- colSums(genes2[ensId$gene_id, ])

qc <- mutate(qc, actb2_f_spike3 = actb2 / f_spike3)

qc %>% ggvis(~group1, ~actb2) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~actb2_f_spike3) %>% ggvis_boxplots()

summary(sort(apply(genes, 2, sum)))
summary(sort(apply(genes2, 2, sum)))

hc1 <- hcluster(t(genes), method = "pearson", link = "average") %>% as.phylo
par(mfrow = c(1, 1)); plot(hc1, edge.width=2, font=2, cex=0.7, label.offset=1e-3, direction="downward")

# data <- data[apply(data, 1, function (x) length(x[x > 10]) > 2), ]
# fc <- sapply(cond, function(idx) rowMeans(genes[, uid == idx])) %>% as.data.frame
# fc <- mutate(fc, WNONP2WIN=WNONP/WIN, MNONP2MIN=MNONP/MIN, WPLM2WIN=WPLM/WIN, MPLM2MIN=MPLM/MIN)

# write.xlsx(as.data.frame(fc), file = "data.xlsx", sheetName = "FOLD_CHANGE", append = T)
# write.table(fc, file = "temp.txt", quote = F)
