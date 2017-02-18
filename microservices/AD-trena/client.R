library(rzmq)
library(jsonlite)
context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://localhost:5557")

i <- 0
while (i < 5) {
   msg = list(cmd="upcase", status="request", callback="handleUpcase", payload=sprintf("hello from R: %d", i))
   msg.json <- toJSON(msg)
   send.raw.string(socket, msg.json)
   s <- receive.string(socket)
   print(fromJSON(s))
   i <- i + 1
   }
