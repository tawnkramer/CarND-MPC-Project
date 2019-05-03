import zmq
import numpy as np
import time
import json
import math

port = 5555
context = zmq.Context()
print ("Connecting to server...")
socket = context.socket(zmq.REQ)
socket.connect ("tcp://localhost:%d" % port)


def make_data(path_len=100):
  d = { 'x' : 0, 'y': -1, 'psi': 0.0, 'speed' : 1.0, 'ptsx' : [], 'ptsy' : [] }

  for i in range(path_len):
      d['ptsx'].append(1.0 * i)
      d['ptsy'].append(math.sin(i * 0.01))

  return json.dumps(d).encode('UTF-8')

state = make_data()
socket.send(state)
message = socket.recv()
obj = json.loads(message.decode('UTF-8'))

i = 0
start = time.time()
numIter = 50
while i < numIter:
  i += 1
  socket.send(state)
  message = socket.recv()
  obj = json.loads(message.decode('UTF-8'))
  print(obj)

dur = time.time() - start
print('%d predictions took %.2f seconds.' % (numIter, dur))
print("Inference at %.2f FPS." % (float(numIter) / dur))
