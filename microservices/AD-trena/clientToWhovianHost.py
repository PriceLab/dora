import zmq
socketContext = zmq.Context()
socket = socketContext.socket(zmq.REQ)
socket.connect("tcp://whovian:%s" % '5556')

for i in range(5):
   socket.send_string("hello from python %d" % i)
   print(socket.recv_string())

