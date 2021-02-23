import sys
import numpy as np
import libsofaTerminal as sofa
import scipy.sparse as sp
import time



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

def doTimer(R):
	NPOINTS=setRes("DeformArmadillo1", R)
	A=time.perf_counter()
	for i in range(1, 100):
		sofa.step()
	B=time.perf_counter()
	return [NPOINTS, B-A]

#####SOFA
sofa.open("UnbuiltMatrix.scn")
[P,T] = doTimer(0.8)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.0)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.2)
print("TIME INCOMING NPTS=", P, " TIME =", T)

sofa.open("IncomingSparseMatrix.scn")
[P,T] = doTimer(0.8)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.0)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.2)
print("TIME INCOMING NPTS=", P, " TIME =", T)

sofa.open("FullIncomingSparseMatrix.scn")
[P,T] = doTimer(0.8)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.0)
print("TIME INCOMING NPTS=", P, " TIME =", T)
[P,T] = doTimer(1.2)
print("TIME INCOMING NPTS=", P, " TIME =", T)


