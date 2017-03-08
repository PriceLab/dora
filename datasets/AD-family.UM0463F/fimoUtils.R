# ~/github/dora/notebooks/AD-family.UM0463F/src/fimoUtils.R
#------------------------------------------------------------------------------------------------------------------------
library(RPostgreSQL)
library(GenomicRanges)
library(TReNA)
library(FimoClient)
library(RUnit)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
dbSNP <- SNPlocs.Hsapiens.dbSNP144.GRCh38
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.snps")){
   tbl.snps <- read.table("~/github/dora/datasets/AD-family.UM0463F/snps59.hg38.bed", sep="\t", as.is=TRUE)
   colnames(tbl.snps) <- c("chrom", "start", "end", "name");
   tokens <- strsplit(tbl.snps$name, "-", fixed=TRUE)
   tbl.snps$wt <- unlist(lapply(tokens, "[", 2))
   tbl.snps$mut <- unlist(lapply(tokens, "[", 3))
   tbl.snps$rsid <- unlist(lapply(tokens, "[", 1))
   }

if(!exists("fimo.service")){
   fimo.service <-  FimoClient("whovian", 5558, quiet=TRUE)
   checkEquals(nrow(requestMatch(fimo.service, list(bogus='NNNNN'))), 0)
   }

if(!exists("chrom.lengths")){
   library(TxDb.Hsapiens.UCSC.hg19.knownGene)   #  version 3.2.2
   tbl.seqInfo <- seqlevels(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
   seqlengths(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))[paste("chr", c(1:22, "X", "Y"), sep="")]
   chrom.lengths <- as.list(seqlengths(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))[paste("chr", c(1:22), sep="")])
  }

db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
tmp <- checkEquals(dbListTables(db.gtf), "hg38human")
if(!exists("tbl.geneInfo")){
   query <- "select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding'"
   tbl.geneInfo <- dbGetQuery(db.gtf, query)
   checkTrue(nrow(tbl.geneInfo) > 19000)
   checkTrue(ncol(tbl.geneInfo) > 25)
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.createWtAndMutSequences()
   test.doComparativeFimo()
   test.toBed9()
   test.intersectWithFootprints()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
doComparativeFimo <- function(chrom, base, wt, mut, snpName, flank, quiet=TRUE)
{
   #wt.sequence <- getSequenceByLoc(dna.service, chrom, base-flank, base+flank)

   ref.sequence <- as.character(getSeq(hg38, chrom, base-flank, base+flank))
   wt.sequence <- paste(substr(ref.sequence, 1, flank), wt, substr(ref.sequence, flank+2, 1+(2*flank)), sep="")
   mut.sequence <- paste(substr(ref.sequence, 1, flank), mut, substr(ref.sequence, flank+2, 1+(2*flank)), sep="")

   if(!quiet){
      printf("--- doComparativeFimo %s:%d %s/%s", chrom, base, wt, mut)
      printf("wt:             %s", wt)
      printf("mut:            %s", mut)
      printf("retrieved:      %s  %s-%s-%s", wt.sequence, substr(wt.sequence, 1, flank),
             substr(wt.sequence, flank+1, flank+1), substr(wt.sequence, flank+2, 1+(2*flank)))
      printf("left flank:     %s", substr(wt.sequence, 1, flank))
      printf("right flank:    %s", substr(wt.sequence, flank+1, nchar(wt.sequence)))
      printf("wt as expected: %s", substr(wt.sequence, flank+1, flank+1) == wt)
      }

   if(mut.sequence == wt.sequence) { # due to the flakey data in tbl.snps.level_1.RData
      printf("    suspicious igap report at %s:%d - mutation same as reference", chrom, base)
      result <- data.frame()
      }
   else{
     query <- list(wt=wt.sequence, mut=mut.sequence)
     result <- requestMatch(fimo.service, query)
     result$X.pattern.name <- as.character(result$X.pattern.name)
     if(!quiet) print(result)
     }

   status <- "noMotif"

   if(nrow(result) == 0){
      if(!quiet) printf("    no motifs in wt or mut")
      }
   else{
     wt.motifs <- unique(subset(result, sequence.name=="wt")$X.pattern.name)
     mut.motifs <- unique(subset(result, sequence.name=="mut")$X.pattern.name)
     novel.mut.motifs <- setdiff(mut.motifs, wt.motifs)
     lost.wt.motifs <- setdiff(wt.motifs, mut.motifs)
     if(!quiet){
       printf("    lost.wt.motifs:   %s", paste(lost.wt.motifs, collapse=","))
       printf("    novel.mut.motifs: %s", paste(novel.mut.motifs, collapse=","))
       }
     if((length(lost.wt.motifs) == 0) & length(novel.mut.motifs) == 0)
         status <- "noChange"
     if((length(lost.wt.motifs) > 0) & length(novel.mut.motifs) > 0)
         status <- "lossAndGain"
     if((length(lost.wt.motifs) > 0) & length(novel.mut.motifs) == 0)
         status <- "loss"
     if((length(lost.wt.motifs) == 0) & length(novel.mut.motifs) > 0)
         status <- "gain"
     } # else: nrow > 0

    if(nrow(result) > 0)
       result$snpName <- snpName

    return(list(chrom=chrom, base=base, wt=wt, mut=mut, flank=flank, status=status, table=result))

} # doComparativeFimo
#------------------------------------------------------------------------------------------------------------------------
test.doComparativeFimo <- function(chrom, base, wt, mut)
{
   printf("--- test.doComparativeFimo")

   x <- with(tbl.snps[4,], doComparativeFimo(chrom, start, wt, mut, rsid, 10))
   checkEquals(sort(names(x)), c("base", "chrom", "flank", "mut", "status", "table", "wt"))
   checkEquals(x$chrom, "chr1")
   checkEquals(x$base, 167474522)

   checkEquals(x$flank, 10)
   checkEquals(x$mut, "G")
   checkEquals(x$wt,  "A")
   checkTrue(is(x$table, "data.frame"))
   checkEquals(x$status, "gain")

} # test.doComparativeFimo
#------------------------------------------------------------------------------------------------------------------------
createWtAndMutSequences <- function(chrom, base, wt.base, mut.base, flank=7)
{
   #wt.sequence <- getSequenceByLoc(dna.service, chrom, base-flank, base+flank)
   wt.sequence <- as.character(getSeq(hg38, chrom, base-flank, base+flank))

   mut.sequence <- paste(substr(wt.sequence, 1, flank), mut.base, substr(wt.sequence, flank+2, 1+(2*flank)), sep="")

   retrieved.wtBase  <- substr(wt.sequence,  flank + 1, flank + 1)
   retrieved.mutBase <- substr(mut.sequence, flank + 1, flank + 1)
   list(wt=wt.sequence, mut=mut.sequence, wtBase=retrieved.wtBase, mutBase=retrieved.mutBase)

} # createWtAndMutSequences
#------------------------------------------------------------------------------------------------------------------------
test.createWtAndMutSequences <- function(chrom, base, flank=7)
{
   printf("--- test.createWtAndMutSequences")
   snp <- tbl.snps[30,]
   chrom <- snp$chrom
   base  <- snp$start
   mut  <-  snp$mut
   wt <- snp$wt
   seqs <- createWtAndMutSequences(chrom, base, wt, mut, flank=3)

   #                         |
   checkEquals(seqs$wt,  "TCAGGCC")
   checkEquals(seqs$mut, "TCAAGCC")
   checkEquals(seqs$wtBase, wt)
   checkEquals(seqs$mutBase, mut)

   seqs <- createWtAndMutSequences(chrom, base, wt, mut, flank=7)
   #                             |
   checkEquals(seqs$wt,  "CAAATCAGGCCTCTG")
   checkEquals(seqs$mut, "CAAATCAAGCCTCTG")
   checkEquals(seqs$wtBase, wt)
   checkEquals(seqs$mutBase, mut)

   for (r in 1:30){
      #printf(" ------ tbl.snps, row %d", r)
      wt  <- tbl.snps$wt[r]
      mut <- tbl.snps$mut[r]
      seqs <- createWtAndMutSequences(tbl.snps$chrom[r], tbl.snps$start[r], wt, mut, flank=3)
      #checkEquals(seqs$wtBase, wt)
      checkEquals(seqs$mutBase, mut)
      } # for r

} # test.createWtAndMutSequences
#------------------------------------------------------------------------------------------------------------------------
toBed9 <- function(tbl)
{
   chrom <- tbl$chrom
   start <- tbl$start - 1
   end   <- tbl$start
   name   <- paste(tbl$rsid, tbl$status, sep="_")
   score <- 1
   strand <- rep(".", nrow(tbl))
   thickStart <- start
   thickEnd   <- end

   colors <- list(noMotif="220,220,220",   # gray
                     loss="255,0,0",       # red
                     gain="0,220,0",       # green
              lossAndGain="255,153,0",     # orange
                 noChange="255,0,255")

   color <- unlist(lapply(tbl$status, function(status) colors[[status]]))
   tbl.out <- data.frame(chrom=chrom, start=start, end=end, name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd, color=color, stringsAsFactors=FALSE)
   tbl.out

} # toBed9
#------------------------------------------------------------------------------------------------------------------------
test.toBed9 <- function()
{
  printf("--- test.toBed9")

  tbl.test <- tbl.snps[1:5,]
  status <- vector(mode="character", length=nrow(tbl.test))

   for(r in 1:nrow(tbl.test)){
      x <- doComparativeFimo(tbl.test$chrom[r], tbl.test$start[r], tbl.test$wt[r],
                             tbl.test$mut[r], tbl.test$rsid[r], 10)
      status[[r]] <- x$status
      }
  tbl.test$status <- unlist(status)


  tbl.b9 <- toBed9(tbl.test)
  checkEquals(dim(tbl.b9), c(nrow(tbl.test), 9))
  checkEquals(colnames(tbl.b9), c("chrom","start","end","name","score","strand","thickStart","thickEnd","color"))
     # spot check first row, 3 values
  with(tbl.b9[1,], checkEquals(name, "rs145175987_loss"),
                   checkEquals(color, "255,0,0"),
                   checkEquals(start, 167461077))

} # test.toBed9
#------------------------------------------------------------------------------------------------------------------------
assessAllSNPs.write.bed9 <- function(outputFile="UM0463_snps_assayed.bed")
{
   status <- vector(mode="character", length=nrow(tbl.snps))
   for(r in 1:nrow(tbl.snps)){
      printf("comparative fimo on snp %d", r)
      x <- doComparativeFimo(tbl.snps$chrom[r], tbl.snps$start[r], tbl.snps$wt[r],
                             tbl.snps$mut[r], tbl.snps$rsid[r], 10)
      status[[r]] <- x$status
      }


  tbl.snps$status <- unlist(status)
  tbl.b9 <- toBed9(tbl.snps)
  checkEquals(dim(tbl.b9), c(nrow(tbl.snps), 9))
  checkEquals(colnames(tbl.b9), c("chrom","start","end","name","score","strand","thickStart","thickEnd","color"))
  printf("writing sample file b9.bed (%d, %d)", nrow(tbl.b9), ncol(tbl.b9))
  write.table(tbl.b9, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, file=outputFile)
  # system("scp b9.bed pshannon@whovian:/local/httpd/vhosts/pshannon/annotations/UM0463_snps_assayed.bed")


} # run.all
#------------------------------------------------------------------------------------------------------------------------
intersectWithFootprints <- function(fp.bed4.file, assayedSnps.bed9.file, outputFile=NULL, shoulder=0)
{
   tbl.fp <- read.table(fp.bed4.file, sep="\t", header=FALSE, as.is=TRUE)
   colnames(tbl.fp) <- c("chrom", "start", "end", "score")
   tbl.snp <- read.table(assayedSnps.bed9.file, sep="\t", header=FALSE, as.is=TRUE)
   colnames(tbl.snp) <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "color")
   gr1 <- with(tbl.snp, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
   gr2 <- with(tbl.fp, GRanges(seqnames=chrom, IRanges(start=start-shoulder, end=end+shoulder)))
   tbl.overlaps <- as.data.frame(findOverlaps(gr1, gr2, type='any'))
   indices <- unique(tbl.overlaps$queryHits)
   printf("number of snps in fps: %d", length(indices))
   tbl.out <- tbl.snp[indices,]

   if(!is.null(outputFile))
      write.table(tbl.out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file=outputFile)
   else
      return(tbl.out)

} # intersectWithFootprints
#------------------------------------------------------------------------------------------------------------------------
test.intersectWithFootprints <- function()
{
   printf("--- test.intersectWithFootprints")
   assayedSnps.bed9.file <- "UM0463_snps_assayed.bed"
   checkTrue(file.exists(assayedSnps.bed9.file))

   fp.bed4.file <- "lizBlue.brain.bed";
   checkTrue(file.exists(fp.bed4.file))
   shoulder <- 0;
   tbl.snpI.brain.0 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=NULL, shoulder=shoulder)
   checkEquals(dim(tbl.snpI.brain.0), c(3, 9))
   checkEquals(tbl.snpI.brain.0$name, c("rs192547848_loss", "rs10918708_gain", "rs75653497_loss"))


   fp.bed4.file <- "lizBlue.lymphoblast.bed";
   checkTrue(file.exists(fp.bed4.file))
   shoulder <- 0;
   tbl.snpI.lymphoblast.0 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=NULL, shoulder=shoulder)
   checkEquals(dim(tbl.snpI.lymphoblast.0), c(4, 9))
   checkEquals(tbl.snpI.lymphoblast.0$name, c("rs192547848_loss", "rs61497249_noMotif", "rs75653497_loss", "rs10918754_noMotif"))


   fp.bed4.file <- "lizBlue.brain.bed";
   checkTrue(file.exists(fp.bed4.file))
   shoulder <- 10;
   tbl.snpI.brain.10 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=NULL, shoulder=shoulder)
   checkEquals(dim(tbl.snpI.brain.10), c(8, 9))
   checkEquals(tbl.snpI.brain.10$name, c("rs114187767_noMotif", "rs192547848_loss", "rs115514202_noMotif", "rs10918708_gain",
                                         "rs61497249_noMotif", "rs17477053_gain", "rs75653497_loss", "rs12745388_noMotif"))

   fp.bed4.file <- "lizBlue.lymphoblast.bed";
   checkTrue(file.exists(fp.bed4.file))
   shoulder <- 10;
   tbl.snpI.lymphoblast.10 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=NULL, shoulder=shoulder)
   checkEquals(dim(tbl.snpI.lymphoblast.10), c(7, 9))
   checkEquals(tbl.snpI.lymphoblast.10$name, c("rs114187767_noMotif", "rs192547848_loss", "rs115514202_noMotif",
                                               "rs79796415_noMotif", "rs61497249_noMotif", "rs75653497_loss", "rs10918754_noMotif"))


} # test.intersectWithFootprints
#------------------------------------------------------------------------------------------------------------------------
createAllAssayedSnpBed9Files <- function(copyToWhovian=FALSE)
{
    filenames.in  <- list(brain.fp.bed="lizBlue.brain.bed",
                          lymphoblast.fp.bed="lizBlue.lymphoblast.bed")
    checkTrue(all(file.exists(as.character(filenames.in))))

    filenames.out <- list(all="UM0463_snps_assayed.bed",
                          brain0="UM0463_brain_snps_0.bed",
                          lymphoblast0="UM0463_lymphoblast_snps_0.bed",
                          brain12="UM0463_brain_snps_12.bed",
                          lymphoblast12="UM0463_lymphoblast_snps_12.bed"
                          )

   outfile <-filenames$all
   assessAllSNPs.write.bed9(outputFile=outfile);
   checkTrue(file.exists(outfile))

   fp.bed4.file <- filenames.in$brain.fp.bed
   shoulder <- 0;
   outfile <- filenames.out$brain0
   tbl.snpI.brain.0 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=outFile, shoulder=shoulder)
   checkTrue(file.exists(outfile))

   fp.bed4.file <- filenames.in$lymphoblast.fp.bed
   shoulder <- 0;
   outfile <- filenames.out$lymphoblast0
   tbl.snpI.lymphoblast.0 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=outfile, shoulder=shoulder)
   checkTrue(file.exists(outfile))

   fp.bed4.file <- filenames.in$brain.fp.bed
   shoulder <- 12;
   outfile <- filenames.out$brain12
   tbl.snpI.brain.12 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=outfile, shoulder=shoulder)
   checkTrue(file.exists(outfile))

   fp.bed4.file <- filenames.in$lymphoblast.fp.bed
   shoulder <- 12;
   outfile <- filenames.out$lymphoblast12
   tbl.snpI.lymphoblast.12 <- intersectWithFootprints(fp.bed4.file, assayedSnps.bed9.file, outputFile=outfile, shoulder=shoulder)
   checkTrue(file.exists(outfile))

    if(copyToWhovian){
       for(f in as.character(filenames.out)){
          cmd <- sprintf("scp %s pshannon@whovian:/local/httpd/vhosts/pshannon/annotations/", f)
          printf("-- about to execute: %s", cmd)
          system(cmd)
          } # for f
       } # if copyToWhovian

} # createAllAssayedSnpBed9Files
#------------------------------------------------------------------------------------------------------------------------
#if(!interactive())
#   runTests()
