import zmq
import json
socketContext = zmq.Context()
socket = socketContext.socket(zmq.REQ)
socket.connect("tcp://localhost:%s" % '5558')

#---------------------------------------------------------------------------
# the most basic test: send 'ping', expect 'pong'
#---------------------------------------------------------------------------
msg = {"cmd": "ping", "status": "request", "callback": "", "payload": ""}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['payload'][0] == 'pong')

#---------------------------------------------------------------------------
# another basic test: send 'upcase' with payload'someLowerCaseWord',
# expect 'SOMELOWERCASEWORD'
#---------------------------------------------------------------------------
msg = {"cmd": "upcase", "status": "request", "callback": "", "payload": "someLowerCaseWord"}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['payload'][0] == 'SOMELOWERCASEWORD')

#---------------------------------------------------------------------------
# the original core command: createGeneModel, with payload targetGene
# and footprintRegion.  send in an intentionally bogus targetGene name
#---------------------------------------------------------------------------
msg = {"cmd": "createGeneModel", "status": "request", "callback": "",
       "payload": {"targetGene": "VEGFabcdbogus", "footprintRegion": "7:101,165,593-101,165,630"}}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['status'][0] == 'error')
assert(response['payload'][0] == 'no expression data for VEGFabcdbogus')

#---------------------------------------------------------------------------
# the original core command: createGeneModel, with payload targetGene
# and footprintRegion.  use a good gene name
#---------------------------------------------------------------------------
msg = {"cmd": "createGeneModel", "status": "request", "callback": "",
       "payload": {"targetGene": "VGF", "footprintRegion": "7:101,165,593-101,165,630"}}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
status = response['status']
payload = response['payload']
assert(len(payload) == 2)
network = payload['network']
footprints = payload['footprints']
assert(network[0][:36] == '{"elements": [ {"data": {"id": "VGF"')
assert(len(footprints[0]) == 3)


