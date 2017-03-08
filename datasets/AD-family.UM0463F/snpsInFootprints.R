# ~/github/dora/datasets/AD-family.UM0463F/snpsInFootprints.R
#------------------------------------------------------------------------------------------------------------------------
library(GenomicRanges)
#------------------------------------------------------------------------------------------------------------------------
tbl.fp <- read.table("lizBlue.brain.bed", sep="\t", header=FALSE, as.is=TRUE)
colnames(tbl.fp) <- c("chrom", "start", "end", "score")
tbl.snp <- read.table("UM0463_snps.hg38.bed", sep="\t", header=FALSE, as.is=TRUE)
colnames(tbl.snp) <- c("chrom", "start", "end")
gr1 <- with(tbl.snp, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
shoulder <- 10
gr2 <- with(tbl.fp, GRanges(seqnames=chrom, IRanges(start=start-shoulder, end=end+shoulder)))
tbl.overlaps <- as.data.frame(findOverlaps(gr1, gr2, type='any'))
indices <- unique(tbl.overlaps$queryHits)
printf("number of snps in fps: %d", length(indices))
write.table(tbl.snp[indices,], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file="UM0463_snps_in_brainFP.bed")
system("scp UM0463_snps_in_brainFP.bed pshannon@whovian:/local/httpd/vhosts/pshannon/annotations/um043f_snps_in_brainFP.bed")
