# createIgapFimoTrack.R
#------------------------------------------------------------------------------------------------------------------------
library(FimoClient)
#library(getDNAClient)
library(RUnit)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38

#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.gwas"))
   load("../tbl.gwas.level_1.RData")
tbl <- tbl.gwas
tbl$score <- -log10(tbl$P)
#if(!exists("dna.service"))
#    dna.service <- getDNAClient("hg38", quiet=TRUE)

if(!exists("fimo.service"))
   fimo.service <-  FimoClient("whovian", 5558, quiet=TRUE)

if(!exists("chrom.lengths")){
   library(TxDb.Hsapiens.UCSC.hg19.knownGene)   #  version 3.2.2
   tbl.seqInfo <- seqlevels(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
   seqlengths(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))[paste("chr", c(1:22, "X", "Y"), sep="")]
   chrom.lengths <- as.list(seqlengths(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))[paste("chr", c(1:22), sep="")])
   }


#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.createWtAndMutSequences()
   test.doComparativeFimo()
   test.toBed9()
   test.geneCentric.bed9.snps()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
doComparativeFimo <- function(chrom, base, wt, mut, flank, quiet=TRUE)
{
   #wt.sequence <- getSequenceByLoc(dna.service, chrom, base-flank, base+flank)
   wt.sequence <- as.character(getSeq(hg38, chrom, base-flank, base+flank))
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

   mut.sequence <- paste(substr(wt.sequence, 1, flank), mut, substr(wt.sequence, flank+2, 1+(2*flank)), sep="")

   if(mut.sequence == wt.sequence) { # due to the flakey data in tbl.gwas.level_1.RData
      printf("    suspicious igap report at %s:%d - mutation same as reference", chrom, base)
      result <- data.frame()
      }
   else{
     query <- list(wt=wt.sequence, mut=mut.sequence)
     result <- requestMatch(fimo.service, query)
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

    status

} # doComparativeFimo
#------------------------------------------------------------------------------------------------------------------------
test.doComparativeFimo <- function(chrom, base, wt, mut)
{
   printf("--- test.doComparativeFimo")
   wt   <- with(tbl[1,], Non_Effect_allele)
   mut  <- with(tbl[1,], Effect_allele)
   chrom <- with(tbl[1,], CHR)
   base <- with(tbl[1,], BP)

   # checkEquals("gain",        doComparativeFimo("chr5", 88469841, "A", "T", 10, quiet=FALSE))
   checkEquals("gain",        doComparativeFimo("chr1", 2312079, "C", "T", 7))
   checkEquals("noMotif",     doComparativeFimo("chr1", 2428478, "T", "C", 7))
   checkEquals("loss",        doComparativeFimo("chr1", 2316546, "T", "C", 7))
   checkEquals("lossAndGain", doComparativeFimo("chr1", 2314000, "T", "C", 7))
   checkEquals("noChange",    doComparativeFimo("chr1", 1068249, "C", "T", 10))

     # some hand-picked examples to bring out all possible results: noMotif, noChange, loss, gain, lossAndGain
   for(r in sample(1:100, size=10)){
     chrom <- tbl$CHR[r]
     base  <- tbl$BP[r]
     mut  <- tbl$Non_Effect_allele[r];
     wt <- tbl$Effect_allele[r];
     printf("%d) %s:%d", r, chrom, base)
     status <- doComparativeFimo(chrom, base, wt, mut, flank=10)
     #printf("status: %s", status)
     #if(nrow(tbl.comp) > 0) print(tbl.comp)
     } # for r

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
   snp <- tbl[30,]
   chrom <- snp$CHR
   base  <- snp$BP
   mut  <-  snp$Effect_allele
   wt <- snp$Non_Effect_allele
   seqs <- createWtAndMutSequences(chrom, base, wt, mut, flank=3)

   #                         |
   checkEquals(seqs$wt,  "GGGCTGG")
   checkEquals(seqs$mut, "GGGTTGG")
   checkEquals(seqs$wtBase, wt)
   checkEquals(seqs$mutBase, mut)

   seqs <- createWtAndMutSequences(chrom, base, wt, mut, flank=7)
   #                             |
   checkEquals(seqs$wt,  "CCTCGGGCTGGGGAC")
   checkEquals(seqs$mut, "CCTCGGGTTGGGGAC")
   checkEquals(seqs$wtBase, wt)
   checkEquals(seqs$mutBase, mut)

   for (r in 1:30){
      #printf(" ------ tbl.gwas, row %d", r)
      wt  <- tbl$Non_Effect_allele[r]
      mut <- tbl$Effect_allele[r]
      seqs <- createWtAndMutSequences(tbl$CHR[r], tbl$BP[r], wt, mut, flank=3)
      #checkEquals(seqs$wtBase, wt)
      checkEquals(seqs$mutBase, mut)
      } # for r

} # test.createWtAndMutSequences
#------------------------------------------------------------------------------------------------------------------------
test.gwasReferenceBases <- function()
{
   printf("--- test.gwasReferenceBases")
   complement <- function(base){
     if(base == "A") return ("T");
     if(base == "T") return ("A");
     if(base == "G") return ("C");
     if(base == "C") return ("G");
     stop("unrecognized base");
     }

   for(r in 1:10){
      claimed.base <- tbl$Non_Effect_allele[r]
      name <- tbl$SNP[r]
      chrom <- tbl$CHR[r]
      loc <- tbl$BP[r]
      #actual.base <- getSequenceByLoc(dna.service, chrom, loc, loc)
      actual.base <-as.character(getSeq(hg38, chrom, loc, loc))
      checkTrue(actual.base %in% c(claimed.base, complement(claimed.base)))
      printf("%s:%d %20s: %s %s %s", chrom, loc, name, claimed.base, actual.base, claimed.base == actual.base)
      } # for r

} # test.gwasReferenceBases
#------------------------------------------------------------------------------------------------------------------------
old.run.mef2c <- function()
{
  chrom <- "chr5"
  start <- 88869312
  end <-   start + 100000
  start <- 88469841 - 20
  end <- 89397039 + 20

  tbl.mef2c <- subset(tbl, CHR==chrom & BP > start & BP < end)   # 67 rows
  motif.status <- vector(mode="character", length=nrow(tbl.mef2c))

  for(r in 1:nrow(tbl.mef2c)){
     chrom <- tbl.mef2c$CHR[r]
     base  <- tbl.mef2c$BP[r]
     mut  <- tbl.mef2c$Non_Effect_allele[r];
     wt <- tbl.mef2c$Effect_allele[r];
     motif.status[r] <- doComparativeFimo(chrom, base, wt, mut, flank=10)
     #printf("status: %s", status)
     #if(nrow(tbl.comp) > 0) print(tbl.comp)
     } # for r

  tbl.mef2c$status <- as.character(motif.status)

  tbl.cory <- read.table("mef2c.corys.snps.tsv", sep="\t", header=FALSE, as.is=TRUE)
  colnames(tbl.cory) <- c("chrom", "start", "end", "rsid", "score", "motif", "ref", "mut")
  motif.status <- vector(mode="character", length=nrow(tbl.cory))

  for(r in 1:nrow(tbl.cory)){
     chrom <- tbl.cory$chrom[r];
     base <- tbl.cory$start[r];
     wt <- tbl.cory$ref[r];
     mut <- tbl.cory$mut[r];
     motif.status[r] <- doComparativeFimo(chrom, base, wt, mut, flank=10)
     }

  tbl.cory$status <- as.character(motif.status)
  list(igap=tbl.mef2c, cory=tbl.cory)

} # old.run.mef2c
#------------------------------------------------------------------------------------------------------------------------
run.chromosome <- function(chrom)
{
  start <- 1
  end <- chrom.lengths[[chrom]]

  tbl.chrom <- subset(tbl, CHR==chrom & BP > start & BP < end)   # 67 rows
  motif.status <- vector(mode="character", length=nrow(tbl.chrom))
  printf("evaluating %d igap snps", nrow(tbl.chrom))

  for(r in 1:nrow(tbl.chrom)){
     chrom <- tbl.chrom$CHR[r]
     base  <- tbl.chrom$BP[r]
     mut  <- tbl.chrom$Non_Effect_allele[r];
     wt <- tbl.chrom$Effect_allele[r];
     motif.status[r] <- doComparativeFimo(chrom, base, wt, mut, flank=10)
     } # for r

  tbl.chrom$status <- as.character(motif.status)

  tbl.b9 <- toBed9(tbl.chrom)
  filename <- sprintf("%ssnps.bed", chrom)
  write.table(tbl.b9, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, file=filename)
  cmd <- sprintf("scp %s pshannon@whovian:/local/httpd/vhosts/pshannon/annotations/", filename)
  system(cmd)

} # run.chromosome
#------------------------------------------------------------------------------------------------------------------------
createBedTracks <- function()
{
  if(!exists("statusPredictions"))
     statusPredictions <<- old.run.mef2c()
  tbl.full <- statusPredictions$igap

  createFile <- function(tbl, selectedStatus){
     tbl.status <- tbl.full[tbl.full$status==selectedStatus,]
     browser()
     filename <- sprintf("snp_%s.bed", selectedStatus)
     write.table(tbl.status[, c(1, 3, 3)], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file=filename);
     }

  createFile(tbl.full, "gain")
  createFile(tbl.full, "loss")
  createFile(tbl.full, "lossAndGain")


} # createBedTracks
#------------------------------------------------------------------------------------------------------------------------
toBed9 <- function(tbl)
{
   chrom <- tbl$CHR
   start <- tbl$BP - 1
   end   <- tbl$BP
   name  <- tbl$SNP
   score <- tbl$score
   strand <- rep(".", nrow(tbl))
   thickStart <- start
   thickEnd   <- end

   colors <- list(noMotif="220,220,220",
                  loss="255,0,0",
                  gain="0,160,0",
                  noChange="255,0,255",
                  lossAndGain="0,0,255")
   color <- unlist(lapply(tbl$status, function(status) colors[[status]]))
   tbl.out <- data.frame(chrom=chrom, start=start, end=end, name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd, color=color, stringsAsFactors=FALSE)
   tbl.out

} # toBed9
#------------------------------------------------------------------------------------------------------------------------
test.toBed9 <- function()
{
  printf("--- test.toBed9")

  if(!exists("statusPredictions"))
    load("statusPredictions.RData")

  tbl.full <- statusPredictions$igap
  stati <- unique(tbl.full$status)
  tbl.test <- tbl.full[match(stati, tbl.full$status),]
  tbl.b9 <- toBed9(tbl.test)
  checkEquals(dim(tbl.b9), c(nrow(tbl.test), 9))
  checkEquals(colnames(tbl.b9), c("chrom","start","end","name","score","strand","thickStart","thickEnd","color"))
  printf("writing sample file b9.bed (%d, %d)", nrow(tbl.b9), ncol(tbl.b9))
  write.table(tbl.b9, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, file="b9.bed")

} # test.toBed9
#------------------------------------------------------------------------------------------------------------------------
explore.trem2 <- function()
{
  chrom <- "chr6"
  trem2.region.start <- 41150938
  trem2.region.end   <- 41196259    # 45k
  tbl.trem2.region <- subset(tbl, BP >= trem2.region.start & BP <= trem2.region.end)   # 102
    # unique rsids provided for each

  max <- nrow(tbl.trem2.region)
  #max <- 10
  motif.status <- vector(mode="character", length=max)
  for(r in 1:max){
     rsid <- tbl.trem2.region$SNP[r]
     loc <- tbl.trem2.region$BP[r]
     ambiguity.code <- snpsById(snps, rsid)$alleles_as_ambig
     elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
     elements <- strsplit(elements.string,'')[[1]]
     wt <- as.character(getSeq(hg38, chrom, loc, loc))
     mut <- setdiff(elements, wt)
     status <- ""
     for(m in 1:length(mut)){
        mut.1 <- mut[m]
        new.status <- doComparativeFimo(chrom, loc, wt, mut.1, 10)
        if(m == 1)
           status <- new.status
        else
          status <- paste(status, new.status, sep=";")
        printf("%3d) %20s %8d: %s -> %s [%s]", r, rsid, loc, wt, mut.1, status)
       }
     motif.status [r] <- status
     } # for r

   tbl.trem2.region$status <- motif.status
   write.table(tbl.trem2.region, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", file="trem2.region.tsv")
   browser();
   tbl.edited <- read.table("trem2.region.tsv", sep="\t", header=TRUE, as.is=TRUE)
   tbl.b9 <- toBed9(tbl.edited)
   write.table(tbl.b9, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, file="trem2.region.bed")


} # explore.trem2
#------------------------------------------------------------------------------------------------------------------------
geneCentric.bed9.snps <- function(chrom, region.start, region.end, gene, manualEdit=FALSE)
{
  tbl.region <- subset(tbl, CHR==chrom & BP >= region.start & BP <= region.end)   # 102
    # unique rsids provided for each
  max <- nrow(tbl.region)
  #browser()
  #max <- 10
  motif.status <- vector(mode="character", length=max)
  for(r in 1:max){
     printf("assessing motif change for %d/%d", r, max);
     rsid <- tbl.region$SNP[r]
     if(rsid =="rs80174570") rsid <- "rs4614309"
     loc <- tbl.region$BP[r]
     ambiguity.code <- snpsById(snps, rsid)$alleles_as_ambig
     elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
     elements <- strsplit(elements.string,'')[[1]]
     wt <- as.character(getSeq(hg38, chrom, loc, loc))
     mut <- setdiff(elements, wt)
     status <- ""
     for(m in 1:length(mut)){
        mut.1 <- mut[m]
        new.status <- doComparativeFimo(chrom, loc, wt, mut.1, 10)
        if(m == 1)
           status <- new.status
        else
          status <- paste(status, new.status, sep=";")
        printf("%3d) %20s %8d: %s -> %s [%s]", r, rsid, loc, wt, mut.1, status)
       }
     motif.status [r] <- status
     #browser()
     #xyz <- 99
     } # for r

   tbl.region$status <- motif.status
   tsv.filename <- sprintf("%s.region.tsv", gene)
   write.table(tbl.region, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", file=tsv.filename)
   if(manualEdit){
      printf("manually edit %s, resolving any plural snp assay values at single loc", tsv.filename)
      browser();
      tbl.edited <- read.table(tsv.filename, sep="\t", header=TRUE, as.is=TRUE)
      tbl.b9 <- toBed9(tbl.edited)
      }
   else{
     tbl.b9 <- tbl.region
      }
   bed.filename <- sprintf("%s.region.bed", gene)
   write.table(tbl.b9, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, file=bed.filename);
   tbl.b9

} # geneCentric.bed9.snps
#------------------------------------------------------------------------------------------------------------------------
test.geneCentric.bed9.snps <- function()
{
   tbl.b9 <- geneCentric.bed9.snps("chr11", 927007, 927125, "AP2A2")
   checkEquals(dim(tbl.b9), c(3, 11))

} # test.geneCentric.bed9.snps
#------------------------------------------------------------------------------------------------------------------------
