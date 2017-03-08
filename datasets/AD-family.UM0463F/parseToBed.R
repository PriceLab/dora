

file <- "UM0463F_BonfSNPs_vep80.tsv"
stopifnot(file.exists(file))
tbl.vep80 <- read.table(file, sep="\t", header=FALSE, as.is=TRUE)
colnames(tbl.vep80) <- c ("chrom", "start", "mut", "str1", "str2",
                          "nearbyGene", "ensg", "enst", "classification",
                          "x10", "x11", "x12", "x13", "x14", "rsid", "misc")
tokens <- strsplit(tbl.vep80$str1, split=":")
ref <- unlist(lapply(tokens, "[", 3))
tbl.vep80$ref <- ref
tbl.vep80$desc <- with(tbl.vep80, paste(rsid, ref, mut, sep="-"))

tbl.hg19.bed <- unique(tbl.vep80[, c("chrom", "start", "start", "desc", "ref", "mut")])
tbl.hg19.bed$chrom <- paste("chr", tbl.hg19.bed$chrom, sep="")
write.table(tbl.hg19.bed[, 1:4], row.names=FALSE, col.names=FALSE, quote=FALSE, file="snps59.bed", sep="\t")

source("~/github/snpFoot/R/liftoverHg19BedFiles.R")
liftoverBedFile.19.38("snps59.bed")
tbl.hg38 <- read.table("snps59.hg38.bed", sep="\t", as.is=TRUE)
colnames(tbl.hg38) <- c("chrom", "start", "endpos", "name")

source("../../notebooks/AD-trena/src/igapGwasFimoUtils.R")

assay <- list()
for(r in 1:nrow(tbl.hg38)){
  ref <- tbl.hg19.bed$ref[r]
  mut <- tbl.hg19.bed$mut[r]
  name <- tbl.hg38$name[r]
  assay[[r]] <- doComparativeFimo(tbl.hg38$chrom[r], tbl.hg38$start[r], ref, mut, 10)
  printf("%20s: %s", name, assay[[r]]$status)
  }

tbl.hg38$status <- unlist(lapply(assay, function(e) e$status))
gainers <- which(lapply(assay, function(x) x$status) == "gain")
library(RPostgreSQL)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="whovian")
checkTrue(dbExistsTable(db, "motifsgenes"))
tbl.mg <- dbGetQuery(db, "select * from motifsgenes")
db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
query <- "select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding'"
tbl <- dbGetQuery(db.gtf, query);     # 19797 x 30

print(load("../AD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData"))
pos <- 167474522;
goi <- subset(tbl, chr=="chr1" & start > pos -10e6 & endpos < pos + 10e6)$gene_name
goi <- intersect(goi, rownames(mtx))

motif.oi <- unique(assay[[4]]$table$X.pattern.name)
tfoi <- unique(subset(tbl.mg, motif %in% motif.oi)$tf)
goi <- unique(c(tfoi, goi))
mtx.sub <- mtx[goi,]
library(gplots)
rc <- rainbow(nrow(mtx), start=0, end=.3)
cc <- rainbow(ncol(mtx), start=0, end=.3)

hm <- heatmap.2(mtx.sub, trace="none", col=cm.colors(255), scale="column",
                #RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                xlab="specification variables", ylab= "Car Models",
                main="",
                density="density")
