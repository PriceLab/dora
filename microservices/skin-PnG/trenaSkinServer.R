# server.R: explore creation and provision of gene models
#------------------------------------------------------------------------------------------------------------------------
PORT=5550
#------------------------------------------------------------------------------------------------------------------------
library(rzmq)
library(jsonlite)
library(RPostgreSQL)
library(TReNA)
library(stringr)
library(graph)
library(RUnit)
library(RCyjs)
source("../trenaCommon/trenaServer.R")
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx.protectedAndExposed")){

   printf("loading expression matrix: %s", load("~/github/dora/datasets/skin/mtx.protectedAndExposed.RData"))
   mtx.protectedAndExposed <- mtx.pAndE
   print(dim(mtx.protectedAndExposed))

   printf("loading expression matrix: %s", load("~/github/dora/datasets/skin/gtex.fib.RData"))
   mtx.gtexFibroblast <- gtex.fib
   print(dim(mtx.gtexFibroblast))

   printf("loading expression matrix: %s", load("~/github/dora/datasets/skin/gtex.primary.RData"))
   mtx.gtexPrimary <- gtex.primary
   print(dim(mtx.gtexPrimary))

   printf("loading expression matrix: %s", load("~/github/dora/datasets/skin/mayo.temporalCortex.RData"))
   print(dim(mtx.temporalCortex))

   printf("loading expression matrix: %s", load("~/github/dora/datasets/skin/mayo.cerebellum.RData"))
   print(dim(mtx.cerebellum))
   }

genome.db.uri    <- "postgres://whovian/hg38"             # has gtf and motifsgenes tables
footprint.db.uri <- "postgres://whovian/skin_hint"        # has hits and regions tables
if(!exists("fpf"))
   fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=TRUE)
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
                            payload=c("skinProtectedAndExposed", "gtexFibroblast", "gtexPrimary",
                                      "mayoTCX", "mayoCER"))
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
           printf("--- executing createGeneModel, with payload:")
           print(msg$payload)              
           targetGene <- msg$payload$targetGene;
           genomicRegion <- msg$payload$footprintRegion
           printf("--- genomicRegion")
           print(genomicRegion)
           expressionMatrixName <- msg$payload$matrix
           mtx.found <- TRUE
           if(expressionMatrixName == "skinProtectedAndExposed")
              mtx <- mtx.protectedAndExposed
           else if(expressionMatrixName == "gtexFibroblast")
               mtx <- mtx.protectedAndExposed <- mtx.gtexFibroblast
           else if(expressionMatrixName == "gtexPrimary")
               mtx <- mtx.gtexPrimary
           else if(expressionMatrixName == "mayoTCX")
               mtx <- mtx.temporalCortex
           else if(expressionMatrixName == "mayoCER")
               mtx <- mtx.cerebellum
           else
              mtx.found <- FALSE
           if(mtx.found) {
              result <- createGeneModel(mtx, targetGene, genomicRegion)
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
