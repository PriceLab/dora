library(TReNA)  # for FootprintFinder
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
roi.chrom <- "chr1"
roi.start <- 167007000
roi.end   <- 168788000
printf("span: %d", roi.end - roi.start)   # 1.78M
#------------------------------------------------------------------------------------------------------------------------
# psql -U trena --host whovian brain_hint
# psql -U trena --host whovian skin_hint
# psql -U trena --host whovian lymphoblast_hint
genome.db.uri    <- "postgres://whovian/hg38"                           # has gtf and motifsgenes tables
brain.footprint.db.uri <- "postgres://whovian/brain_hint"               # has hits and regions tables
skin.footprint.db.uri <- "postgres://whovian/skin_hint"                 # has hits and regions tables
lymphoblast.footprint.db.uri <- "postgres://whovian/lymphoblast_hint"   # has hits and regions tables

if(!exists("fpf.brain"))
   fpf.brain <- FootprintFinder(genome.db.uri, brain.footprint.db.uri, quiet=FALSE)

if(!exists("fpf.skin"))
   fpf.skin <- FootprintFinder(genome.db.uri, skin.footprint.db.uri, quiet=FALSE)

if(!exists("fpf.lymphoblast"))
   fpf.lymphoblast <- FootprintFinder(genome.db.uri, lymphoblast.footprint.db.uri, quiet=FALSE)

#------------------------------------------------------------------------------------------------------------------------
# reduce the footprint table to a bed format, with just one row per unique loc
# and a score which counts the number of samples at that exact loc.
# the incoming tbl.fp will also have a unique row for every motif assigned
# to the loc, and is sometimes further inflated by palindromic motifs
# mapped to both strands.  ignore all of that.
condense.footprints <- function(tbl.fp)
{
   tbl.fp.uniq <- unique(tbl.fp[, c("loc", "sample_id")])  # the only two relevant fields
   tbl.reps <- as.data.frame(table(tbl.fp.uniq$loc), stringsAsFactors=FALSE)  # find their distrution
   tbl.repsOfInterest <- subset(tbl.reps, Freq > 0)  # get rid of all singleton locs
   indices <- match(tbl.repsOfInterest$Var1, tbl.fp$loc)   # find out their first occurence
   tbl.chromStartEnd <- tbl.fp[indices, c("chrom", "start", "endpos")]   # extract chrom, start, end
   tbl.out <- cbind(tbl.repsOfInterest, tbl.chromStartEnd)     # broaden loc, sampleCount with chrom,start,end
   colnames(tbl.out) <- c("loc", "score", "chrom", "start", "end")   # fix the colnames

   tbl.out[, c("chrom", "start", "end", "score")]

} # condense.footprints
#------------------------------------------------------------------------------------------------------------------------
test.condense.footprints <- function()
{
   printf("--- test.condense.footprints")
   if(!exists("tbl.fp22"))
      load("tbl.fp.chr22.RData", envir=.GlobalEnv())
   tbl.test <- tbl.fp22; # [1:100,]  # has a 2-dup, a 4-dup and 5 singles
   tbl.condensed <- condense.footprints(tbl.test)

   checkTrue(nrow(tbl.condensed) < 80000)
   checkEquals(ncol(tbl.condensed), 4)
   checkEquals(colnames(tbl.condensed), c("chrom", "start", "end", "score"))

} # test.condense.footprints
#------------------------------------------------------------------------------------------------------------------------
saveBed <- function(tbl.fp, filename)
{
   dim(tbl.fp)   # 4105 17
   tbl.bed <- tbl.fp[, c("chrom", "start", "endpos", "name", "score3")]
   write.table(tbl.bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file="test.bed")


} # saveBed
#------------------------------------------------------------------------------------------------------------------------
iterate <- function()
{
   fpfs <- list(brain=fpf.brain, skin=fpf.skin, lymphoblast=fpf.lymphoblast)
   fpfs <- list(lymphoblast=fpf.lymphoblast)

   for(fp.name in names(fpfs)){
     fpf <- fpfs[[fp.name]]
     tbl.fp <- getFootprintsInRegion(fpf, roi.chrom, roi.start, roi.end)
     tbl.fpReduced <- condense.footprints(tbl.fp)
     filename <- sprintf("lizBlue.%s.bed", fp.name)
     write.table(tbl.fpReduced, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file=filename)
     }

} # iterate
#------------------------------------------------------------------------------------------------------------------------

