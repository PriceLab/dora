# 11:923,951-927,102
#------------------------------------------------------------------------------------------------------------------------
library(TReNA)
library(plyr)
library(FimoClient)
library(RUnit)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
dbSNP <- SNPlocs.Hsapiens.dbSNP144.GRCh38


stopifnot(packageVersion("TReNA") >= '0.99.40')

if(!exists("fimo.services"))
   fimo.service <-  FimoClient("whovian", 5558, quiet=TRUE)
result <- requestMatch(fimo.service, list(bogus='xxxxx'))
tmp <- checkEquals(result, data.frame())

# create tbl.geneInfo for TSS
library(RPostgreSQL)
if(!exists("db.gtf")){
   db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
   query <- "select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding'"
   tbl.geneInfo <- dbGetQuery(db.gtf, query)[, c("chr", "gene_name", "start", "strand")]
   }

if(!exists("mtx.expression")){
  print(load("~/github/dora/datasets/AD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData"))
  mtx.expression <- asinh(mtx)
  }

stopifnot(dim(mtx.expression) == c(25031, 638))

if(!exists("createGeneModel"))
    source("~/github/dora/microservices/AD-trena/server.R")

if(!exists("tbl.gwas")){
  print(load("~/github/dora/datasets/AD/tbl.gwas.level_1.RData"))
  }

if(!exists("doComparativeFimo")){
   source("~/github/dora/notebooks/AD-trena/src/createIgapFimoTrack.R")
   }



#----------------------------------------------------------------------------------------------------
# use the standard microservice TReNA
model.laptm5 <- function()
{
   target.gene <- "LAPTM5"
    # 1:30,753,461-30,759,462
     # subset(tbl.geneInfo, gene_name==target.gene)$start  TODO - it's the end
   tss <- 30757828;
   start <- tss - 4200
   end   <- tss + 1300
   chrom <- "chr1"
   region <- sprintf("%s:%d-%d", chrom, start, end)
   result <- createGeneModel(target.gene, region)
   tbl.gm <- result$tbl

} # model.laptm5
#----------------------------------------------------------------------------------------------------
find.newTFs <- function()
{
     # "rs10794342" "rs10794343" "rs10902233"
  tss <- 30757828;
  start <- tss - 4200
  end   <- tss + 1300
  chrom <- "chr1"

  tbl.prospects <- subset(tbl.gwas, CHR=="chr1" & (BP > start) & (BP < end))

  gene <- "LAPTM5"
  for(r in 1:nrow(tbl.prospects)){
    rsid <- tbl.prospects$SNP[r]
    if(rsid == "rs79037040") rsid <- "rs13278062"
    chrom <- tbl.prospects$CHR[r]
    loc <- tbl.prospects$BP[r]
    #gene <- tbl.prospects$gene_name[r]
    ambiguity.code <- snpsById(dbSNP, rsid)$alleles_as_ambig
    elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
    elements <- strsplit(elements.string,'')[[1]]
    wt <- as.character(getSeq(hg38, chrom, loc, loc))
    mut <- setdiff(elements, wt)
    status <- doComparativeFimo(chrom, loc, wt, mut, 10, quiet=FALSE)
    printf("---- %s, %s: %s", rsid, gene, status)
    }


} # find.newTFs
#----------------------------------------------------------------------------------------------------
# from the !quiet print trace of doComparativeFimo, just above
#
# --------------- doComparativeFimo chr11:924904 C/T
#
# wt:             C
# mut:            T
# retrieved:      CTCTTCCTGGCCGGGACCTGG  CTCTTCCTGG-C-CGGGACCTGG
# left flank:     CTCTTCCTGG
# right flank:    CCGGGACCTGG
# wt as expected: TRUE
#   pattern.name sequence.name start stop strand    score  p.value q.value
# 1       MA0136.2           mut     1   11      - 12.21820 3.11e-05 0.00137
# 2       MA0598.2           mut     1   12      -  9.59016 9.06e-05 0.00362
#
# --- doComparativeFimo chr11:925493 C/T
#
# wt:             C
# mut:            T
# retrieved:      GGCAAGCTAGCGGGGCCGAGG  GGCAAGCTAG-C-GGGGCCGAGG
# left flank:     GGCAAGCTAG
# right flank:    CGGGGCCGAGG
# wt as expected: TRUE
#         X.pattern.name sequence.name start stop strand   score  p.value q.value
# 1 sci09.v2_Plagl1_0972           mut    10   19      + 9.81818 9.47e-05 0.00455
#
# --- doComparativeFimo chr11:925510 A/G
#
# wt:             A
# mut:            G
# retrieved:      GAGGCCCAGGACGAGGCCGTA  GAGGCCCAGG-A-CGAGGCCGTA
# left flank:     GAGGCCCAGG
# right flank:    ACGAGGCCGTA
# wt as expected: TRUE
#    X.pattern.name sequence.name start stop strand    score  p.value  q.value
# 1       ZNF148.p2           mut     1   16      + 13.94740 6.95e-06 0.000167
# 2        MA0146.2            wt     1   14      - 13.60610 1.52e-05 0.000248
# 3        MA0146.2           mut     1   14      - 13.57580 1.55e-05 0.000248
# 4        MA0163.1           mut     1   14      +  7.93878 3.09e-05 0.000988
# 5        MA0810.1           mut     2   13      - 11.92730 4.36e-05 0.001740
# 6        MA0599.1           mut     8   17      - 11.91840 4.63e-05 0.002220
# 7        MA0524.2           mut     2   13      + 11.59180 5.14e-05 0.001390
# 8        MA0811.1           mut     2   13      + 11.61820 5.26e-05 0.001360
# 9        MA0811.1           mut     2   13      - 11.21820 6.81e-05 0.001360
# 10       MA0524.2           mut     2   13      - 11.15310 6.93e-05 0.001390
# 11       MA0039.2           mut     8   17      + 10.44900 7.72e-05 0.003700
#------------------------------------------------------------------------------------------------------------------------
characterize.gained.tfs <- function()
{
    motifs <- c("MA0526.1", "MA0009.2", "MA0147.2", "MA0616.1", "MA0804.1", "MA0649.1")

    if(!exists("tb.mg")){
       db.mg <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="whovian")
       tbl.mg <<- dbGetQuery(db.mg, "select * from motifsgenes")
       }

     new.tfs <- unique(subset(tbl.mg, motif %in% motifs)$tf)  # 62!
     new.known.tfs <- c(intersect(new.tfs, rownames(mtx.expression)), "FOXK1", "HMGA1", "FOXP4")
     cor.list <- lapply(new.known.tfs, function(tf) cor(mtx.expression[tf,], mtx.expression["AP2A2",]))
     names(cor.list) <- new.known.tfs
     tbl.cor <- data.frame(tf=names(cor.list), cor=as.numeric(cor.list), stringsAsFactors=FALSE)
     tbl.cor <- tbl.cor[order(abs(tbl.cor$cor), decreasing=TRUE),]

} # characterize.gained.tfs
#------------------------------------------------------------------------------------------------------------------------
#runTReNA <- function(gene, mtx, candidate.tfs)
#{
#  trena.lasso <- TReNA(mtx, solver="lasso")
#  trena.bs    <- TReNA(mtx, solver="bayesSpike")
#  trena.rf    <- TReNA(mtx, solver="randomForest")
#
#  tfs <- intersect(candidate.tfs, rownames(mtx))
#  printf("%8s: %d tf candidates", gene, length(tfs))
#  tbl.lasso <- solve(trena.lasso, gene, tfs, extraArgs = list(alpha = 1.0))
#  tbl.lasso <- subset(tbl.lasso, abs(beta) > 0.01)
#  tbl.bs <- solve(trena.bs, gene, tfs)
#  if(nrow(tbl.bs) > 0)
#    tbl.bs <- subset(tbl.bs, pval < 0.05)
#
#  suppressWarnings(
#    rf.out <- solve(trena.rf, gene, tfs)
#    )
#
#  tbl.rf = rf.out$edges
#  tbl.rf <-subset(tbl.rf, IncNodePurity >= fivenum(tbl.rf$IncNodePurity)[4])
#  tbl.rf <- tbl.rf[order(tbl.rf$IncNodePurity, decreasing=TRUE),,drop=FALSE]
#
#  tbl.lasso$gene <- rownames(tbl.lasso)
#  tbl.lasso$method <- "lasso"
#  tbl.lasso$score <- tbl.lasso$beta
#
#  tbl.bs$gene <- rownames(tbl.bs)
#  tbl.bs$method <- "bayesSpike"
#  tbl.bs$score <- tbl.bs$beta
#
#  tbl.rf$gene <- rownames(tbl.rf)
#  tbl.rf$method <- "randomForest"
#  tbl.rf$score <- tbl.rf$IncNodePurity
#
#  tbl.all <- rbind.fill(tbl.lasso, tbl.bs, tbl.rf)
#  tbl.all[order(abs(tbl.all$gene.cor), decreasing=TRUE),]
#
#} # runTReNA
##----------------------------------------------------------------------------------------------------
## our demo gene is MEF2C, for which a small rna expression matrix is included
## in the TReNA package:
##
##           gene_id gene_name  chr    start   endpos strand
## 1 ENSG00000081189     MEF2C chr5 88717117 88904257      -
##
goi <- "AP2A2"
tss <- subset(tbl.geneInfo, gene_name==goi)$start  # [1] 924894, on + strand
start <- tss - 1000
end   <- tss + 1000
chrom <- "chr11"

   # just keep rows with some variance across samples
# sd <- apply(mtx.expression, 1, sd)
# deleters <- which(sd < 1)
# if(length(deleters) > 0)
#   mtx.sub <- mtx.sub[-deleters,]
# mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
# mtx.sub.log2 <- log2(mtx.tmp)
# fivenum(mtx.sub.log2)
# genes <- rownames(mtx.sub.log2)
#
#    # no TF self-loops, but target must be in the expression matrix
# stopifnot("MEF2C" %in% genes)
# tf.candidates <- genes[-grep("MEF2C", genes)]
# tbl.all <- runTReNA(goi, mtx.sub.log2, tf.candidates)
# print(head(tbl.all, n=10))
#    # compare the methods
# print(head(table(tbl.all$gene, tbl.all$method)))

