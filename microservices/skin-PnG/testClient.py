# testClient.py:  exercise some of the capabilities offered by the TReNA gene model server
#------------------------------------------------------------------------------------------
PORT = 5558;
import zmq
import json
socketContext = zmq.Context()
socket = socketContext.socket(zmq.REQ)
socket.connect("tcp://localhost:%s" % PORT)
#---------------------------------------------------------------------------
# the most basic test: send 'ping', expect 'pong'
#---------------------------------------------------------------------------
msg = {"cmd": "ping", "status": "request", "callback": "", "payload": ""}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['payload'] == 'pong')

#---------------------------------------------------------------------------
# another basic test: send 'upcase' with payload'someLowerCaseWord',
# expect 'SOMELOWERCASEWORD'
#---------------------------------------------------------------------------
msg = {"cmd": "upcase", "status": "request", "callback": "", "payload": "someLowerCaseWord"}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['payload'] == 'SOMELOWERCASEWORD')

#---------------------------------------------------------------------------
# the TReNA server will often have multiple expression matrices
# to choose from.  return their names here upon request
#---------------------------------------------------------------------------
msg = {"cmd": "getExpressionMatrixNames", "status": "request", "callback": "", "payload": ""}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['status'] == 'success')
payload = response['payload']
assert(response['payload'] == ['skinProtectedAndExposed', 'gtexFibroblast', 'gtexPrimary'])

#---------------------------------------------------------------------------
# the original core command: createGeneModel, with payload targetGene
# and footprintRegion.  send in an intentionally bogus targetGene name
#---------------------------------------------------------------------------
msg = {"cmd": "createGeneModel", "status": "request", "callback": "",
       "payload": {"targetGene": "VEGFabcdbogus",
                   "matrix": "gtexFibroblast",
                   "footprintRegion": "7:101,165,593-101,165,630"}}
msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
assert(response['status'] == 'error')
assert(response['payload'] == 'no expression data for VEGFabcdbogus')

#---------------------------------------------------------------------------
# the original core command: createGeneModel, with payload targetGene
# and footprintRegion.  use a good gene name
#---------------------------------------------------------------------------
msg = {"cmd": "createGeneModel", "status": "request", "callback": "",
       "payload": {"targetGene": "COL1A1",
                   "matrix": "gtexFibroblast",
                   "footprintRegion": "17:50,203,128-50,203,174"}}

msg_json = json.dumps(msg)
socket.send_string(msg_json)
response = json.loads(socket.recv_string())
status = response['status']
payload = response['payload']
assert(list(payload.keys()) == ['network', 'model', 'footprints'])
network = payload['network']
footprints = payload['footprints']
model = payload['model']
assert(network[:40] == '{"elements": [ {"data": {"id": "COL1A1",')
assert(len(footprints[0]) == 17)


