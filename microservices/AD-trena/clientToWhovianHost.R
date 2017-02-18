library(rzmq)
library(jsonlite)
context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://whovian:5556")

i <- 0
while (i < 5) {
   send.raw.string(socket, sprintf("hello from R: %d", i))
   s <- receive.string(socket)
   print(s)
   i <- i + 1
   }
