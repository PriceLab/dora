# trenaBrainServer.R: explore creation and provision of gene models
#------------------------------------------------------------------------------------------------------------------------
source("../trenaCommon/trenaServer.R")
#------------------------------------------------------------------------------------------------------------------------
library(rzmq)
PORT=5551
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx.rosmap")){
   load("../../datasets/AD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
   mtx.asinh <- asinh(mtx)
   medians <- apply(mtx.asinh, 1, median)
   deleters <- as.integer(which(medians <= 0.1))
   mtx.rosmap <- mtx.asinh[-deleters,]
   }
genome.db.uri    <- "postgres://whovian/hg38"             # has gtf and motifsgenes tables
footprint.db.uri <- "postgres://whovian/brain_hint"        # has hits and regions tables
if(!exists("fpf"))
   fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=TRUE)

#------------------------------------------------------------------------------------------------------------------------
trenaBrainServerTests <- function()
{
   test.createGeneModel.brain()

} # trenaBrainServerTests
#------------------------------------------------------------------------------------------------------------------------
demo.rs34423320.mpzl1 <- function()
{
   printf("--- demo.rs34423320")

} # demo.rs34423320
#------------------------------------------------------------------------------------------------------------------------
test.createGeneModel.brain <- function()
{
   printf("--- test.createGeneModel.brain")

          # on the + strand, a large region
   region <- "7:101,165,571-101,166,620"
   target.gene <- "VGF"
   result <- createGeneModel(mtx=mtx.rosmap, target.gene, region)
   tbl.gm <- result$model
   tbl.reg <- result$regulatoryRegions
   msg <- result$msg
   checkEquals(msg, "112 putative TFs found")
   checkTrue(ncol(tbl.gm) >= 10)
   checkTrue(nrow(tbl.gm) >= 8)
   checkTrue("EGR4" %in% tbl.gm$gene)

      # try a gene on the minus strand
   target.gene <- "MEF2C"
   region <- "5:88,904,000-88,909,000"
   result <- createGeneModel(mtx=mtx.rosmap, target.gene, region)
   tbl.gm <- result$model
   checkTrue(ncol(tbl.gm) >= 10)
   checkTrue(nrow(tbl.gm) >= 8)
   checkTrue("PAX7" %in% tbl.gm$gene)

} # test.createGeneModel.brain
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) {

   context = init.context()
   socket = init.socket(context,"ZMQ_REP")
   bind.socket(socket, sprintf("tcp://*:%d", PORT))

   errorFunction <- function(condition){
     printf("==== exception caught ===")
     print(as.character(condition))
     response <- list(cmd=msg$callack, status="error", callback="", payload=as.character(condition));
     send.raw.string(socket, toJSON(response))
     };

   while(TRUE) {
      tryCatch({
        printf("top of receive/send loop")
        raw.message <- receive.string(socket)
        msg = fromJSON(raw.message)
        printf("cmd: %s", msg$cmd)
        print(msg)
        if(msg$cmd == "ping") {
            response <- list(cmd=msg$callack, status="result", callback="", payload="pong")
            }
        else if(msg$cmd == "upcase") {
            response <- list(cmd=msg$callack, status="result", callback="", payload=toupper(msg$payload))
            }
        else if(msg$cmd == "getTestNetwork"){
           infile <- file("vgfModel.json")
           graphModel <- fromJSON(readLines(infile))
           response <- list(cmd=msg$callback, status="result", callback="", payload=graphModel)
           }
        else if(msg$cmd == "getExpressionMatrixNames"){
           printf("getExpressionMatrixNames");
           response <- list(cmd=msg$callback, status="success", callback="",
                            payload=c("rosmap"))
           }
        else if(msg$cmd == "getFootprintsInRegion"){
          footprintRegion <- msg$payload$footprintRegion;
          region.parsed <- extractChromStartEndFromChromLocString(footprintegion)
          chrom <- region.parsed$chrom
          start <- region.parsed$start
          end <-   region.parsed$end
          printf("region parsed: %s:%d-%d", chrom, start, end);
          tbl.fp <- getFootprintsInRegion(fpf, chrom, start, end)
          commandStatus <- "success"
          if(nrow(tbl.f) == 0)
             commandStatus = "failure"
          response <- list(cmd=msg$callback, status=commandStatus, callback="", payload=tbl.fp)
          }
        else if(msg$cmd == "createGeneModel"){
           print(1)
           targetGene <- msg$payload$targetGene;
           print(2)
           genomicRegions <- msg$payload$genomicRegions
           printf("--- genomicRegions")
           print(genomicRegions)
           print(3)
           expressionMatrixName <- msg$payload$matrix
           mtx.found <- TRUE
           if(expressionMatrixName == "rosmap")
              mtx <- mtx.rosmap
           else
              mtx.found <- FALSE
           if(mtx.found) {
              result <- createGeneModel(mtx, targetGene, genomicRegions)
              tbl.gm <- result$model
              tbl.reg <- result$regulatoryRegions
              message <- result$msg
              if(nrow(tbl.gm) == 0){
                response <- list(cmd=msg$callback, status="error", callback="", payload=message)
                }
              else{
                tbl.list <- list(tbl.gm)
                print(5)
                names(tbl.list) <- targetGene
                g <- tablesToFullGraph(tbl.gm, tbl.reg)
                g.lo <- addGeneModelLayout(g)
                g.json <- graphToJSON(g.lo)
                #g.json <- sprintf("network = %s", g.json)
                print(8)
                payload <- list(network=g.json, model=tbl.gm, footprints=tbl.reg)
                response <- list(cmd=msg$callback, status="success", callback="", payload=payload)
                }
              } # mtx.found
           if(!mtx.found){
               response <- list(cmd=msg$callback, status="error", callback="",
                                payload=sprintf("unrecognized matrix name: '%s'",  expressionMatrixName))
               }
           } # createGeneModel
        else {
           response <- list(cmd="handleUnrecognizedCommand", status="error", callback="", payload=toupper(raw.message))
           }
        printf("--- about to send.raw.string")
        json.string <- toJSON(response, dataframe="values")
        # printf("json.string: %s", json.string)
        send.raw.string(socket, json.string)
        Sys.sleep(1)
        }, error=errorFunction); # tryCatch
     } # while (TRUE)

} # if !interactive()
#------------------------------------------------------------------------------------------------------------------------
