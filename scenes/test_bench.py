import sys
import numpy as np
import libsofaTerminal as sofa
import scipy.sparse as sp
import time
import matplotlib.pyplot as plt


def setRes(NODE_NAME,RES):
	ROOT=sofa.root()
	
	GRID=ROOT.getObject("Topology/Grid")
	GRID.clear()
	GRID.getData("cellWidth").setValue(RES)
	GRID.init()

	MAPPING=ROOT.getObject("Topology/Mapping")
	MAPPING.init()
	
	CONTAINER=ROOT.getObject("Topology/Container")
	NPOINS=CONTAINER.getData("nbPoints").getValue()
	POS=CONTAINER.getData("position").getValue(False)
	
	CTX = ROOT.getContext(NODE_NAME)
	
	STATE=CTX.getState()
	STATE.resize(0)
	STATE.resize(NPOINS)
	STATE.getData("rest_position").setValue(POS)
	STATE.getData("reset_position").setValue(POS)
	STATE.reset()
	
	CTX.propagateTopologicalChanges()
	
	return NPOINS

def AssemblyTimer(R):
	NPOINTS=setRes("DeformArmadillo1", R)

	# for i in range(1, 100):
	# 	sofa.step()
	# print(sofa.timer().getTime("MBKBuild"))
	sofa.step()
	sofa.timer().activate()

	for i in range(1, 101):
		sofa.step()
	Assembly=sofa.timer().getTime("MBKBuild")
	print(Assembly)
	# print(sofa.timer().getVal("PCG iterations"))
	sofa.timer().clear()
	return [NPOINTS, Assembly[5]]


input = [1.2, 1.1, 1.0, 0.9, 0.8]
x=[]
Unbuilt = []
LocalIncoming = []
FullIncoming = []
CudaFullIncoming = []

#####SOFA
for i in input:
	sofa.open("UnbuiltMatrix.scn")
	[P,T] = AssemblyTimer(i)
	Unbuilt.append(T)
	x.append(P)

for i in input:
	sofa.open("IncomingSparseMatrix.scn")
	[P,T] = AssemblyTimer(i)
	LocalIncoming.append(T)

for i in input:
	sofa.open("FullIncomingSparseMatrix.scn")
	[P,T] = AssemblyTimer(i)
	FullIncoming.append(T)

for i in input:
	sofa.open("CudaFullIncomingSparseMatrix.scn")
	[P,T] = AssemblyTimer(i)
	CudaFullIncoming.append(T)

print("x =",x)
print("Unbuilt =",Unbuilt)
print("LocalIncoming =",LocalIncoming)
print("FullIncoming =",FullIncoming)
print("CudaFullIncoming =",CudaFullIncoming)

# plt.plot(x, Unbuilt, label = "Unbuilt")
plt.plot(x, LocalIncoming, label = "LocalIncoming")
plt.plot(x, FullIncoming, label = "FullIncoming")
plt.plot(x, CudaFullIncoming, label = "CudaFullIncoming")

plt.xlabel('NB points')
plt.ylabel('time (ms)')
plt.title('comparison assembly ')
# show a legend on the plot
plt.legend()
# Display a figure.
plt.show()
