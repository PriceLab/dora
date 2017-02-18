library(RUnit)
library(rzmq)
library(jsonlite)
#------------------------------------------------------------------------------------------------------------------------
context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://localhost:5557")
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
    test.ping()
    test.createGeneModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test.ping <- function()
{
   printf("--- test.ping")
   msg = list(cmd="ping", status="request", callback="", payload="");
   msg.json <- toJSON(msg)
   send.raw.string(socket, msg.json)
   response <- fromJSON(receive.string(socket))
   checkEquals(response$payload, "pong")

} # test.ping
#------------------------------------------------------------------------------------------------------------------------
test.createGeneModel <- function()
{
   printf("--- test.createGeneModel")

   msg = list(cmd="createGeneModel", status="request", callback="displayGeneModel",
               payload=list(targetGene="VGF", footprintRegion="7:101,165,577-101,165,615"))
   msg.json <- toJSON(msg)
   send.raw.string(socket, msg.json)
   response <- fromJSON(receive.string(socket))
   browser()
   x <- 99


} # test.createGeneModel
#------------------------------------------------------------------------------------------------------------------------
