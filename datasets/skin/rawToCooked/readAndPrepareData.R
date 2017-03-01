library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!file.exists("gtex.rpkm.RData")){
   gtex.meta <- read.table("/local/Cory/gtex/data/GTEx_Data_V6_Annotations_SampleAttributesDS2.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
   gtex.rpkm <- read.table("/local/Cory/gtex/data/GTEX_rpkm2", header=TRUE, stringsAsFactors=FALSE, sep="\t")

   save(gtex.meta, file="gtex.meta.RData")
   save(gtex.rpkm, file="gtex.rpkm.RData")
   }
else{
   load("gtex.meta.RData")
   load("gtex.rpkm.RData")
   }

#------------------------------------------------------------------------------------------------------------------------
# eliminate duplicates
dups <- gtex.rpkm[duplicated(gtex.rpkm$Description),]
tbl.noDups <- gtex.rpkm[!(duplicated(gtex.rpkm$Description)),]
rownames(tbl.noDups) <- tbl.noDups$Description

   # create the data.frame for just samples with protected (unexposed) skin
protected.key <- "Skin - Not Sun Exposed (Suprapubic)"
protected.cols <- grep(protected.key, gtex.meta$SMTSD, fixed=TRUE)
tbl.protected <- gtex.meta[protected.cols, ]   # 273 64
sampleIDs.protected <- gsub("-", ".", tbl.protected$SAMPID, fixed=TRUE)
rpkm.column.names.for.protected <- intersect(sampleIDs.protected, colnames(tbl.noDups))
tbl.protected <- tbl.noDups[, rpkm.column.names.for.protected]
checkEquals(dim(tbl.protected), c(54301, 250))

exposed.key <- "Skin - Sun Exposed (Lower leg)"
exposed.cols <- grep(exposed.key, gtex.meta$SMTSD, fixed=TRUE) # 468
tbl.exposed <- gtex.meta[exposed.cols, ]   # 468 64
sampleIDs.exposed <- gsub("-", ".", tbl.exposed$SAMPID, fixed=TRUE)
rpkm.column.names.for.exposed <- intersect(sampleIDs.exposed, colnames(tbl.noDups))
length(rpkm.column.names.for.exposed)
tbl.exposed <- tbl.noDups[, rpkm.column.names.for.exposed]
checkEquals(dim(tbl.exposed), c(54301, 357))

checkEquals(length(intersect(colnames(tbl.exposed), colnames(tbl.protected))), 0)
tbl.protectedAndExposed <- cbind(tbl.protected, tbl.exposed)
dim(tbl.protectedAndExposed)  #  54301   607
save(tbl.protectedAndExposed, file="tbl.protectedAndExposed.RData")

# exposed <- gtex.meta[which(gtex.meta$SMTSD=="Skin - Sun Exposed (Lower leg)"),]
# fibro <- gtex.meta[which(gtex.meta$SMTSD=="Cells - Transformed fibroblasts"),]
# all_skin <- rbind(protected, exposed, fibro)
# primary <- rbind(protected, exposed)
# 
# selected.columns <- match(head(sampleIDs.protected), colnames(solo))
# 
# rownames(solo) <- solo$Description
# solo$Description <- NULL
# solo$Name <- NULL
# 
# gtex.prot <- as.matrix(solo[,which(protected$SAMPID %in% colnames(solo))])
# gtex.ex <- as.matrix(solo[,which(exposed$SAMPID %in% colnames(solo))])
# gtex.fib <- as.matrix(solo[,which(fibro$SAMPID %in% colnames(solo))])
# gtex.all_skin <- as.matrix(solo[,which(all_skin$SAMPID %in% colnames(solo))])
# gtex.protectedAndExposed <- as.matrix(solo[,which(primary$SAMPID %in% colnames(solo))])
# save(gtex.prot, gtex.ex, gtex.fib, gtex.all_skin, gtex.primary, file="gtex.matrices.allByCondition.RData")
# 
