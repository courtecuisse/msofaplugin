import sys
import numpy as np
import libsofaTerminal as sofa
import scipy.sparse as sp
import time



def setRes(R):
	ROOT=sofa.root()
	
	GRID=ROOT.getObject("SparseGrid")
	GRID.clear()
	GRID.getData("cellWidth").setValue(R)
	GRID.init()

	MAPPING=ROOT.getObject("DeformArmadillo1/mapping")
	MAPPING.init()
	
	CONTAINER=ROOT.getObject("DeformArmadillo1/Container")
	NPOINS=CONTAINER.getData("nbPoints").getValue()
	POS=CONTAINER.getData("position").getValue(False)
	
	STATE=ROOT.getObject("DeformArmadillo1/mstate")
	STATE.resize(0)
	STATE.resize(NPOINS)
	STATE.getData("rest_position").setValue(POS)
	STATE.getData("reset_position").setValue(POS)
	STATE.reset()

	MASS=ROOT.getObject("DeformArmadillo1/mass")
	MASS.getData("indices").setValue(CONTAINER.getData("points").getValue())
	
	ROOT.getContext("DeformArmadillo1").propagateTopologicalChanges()
	
	return NPOINS

def doTimer(R):
	NPOINTS=setRes(R)
	A=time.perf_counter()
	for i in range(1, 100):
		sofa.step()
	B=time.perf_counter()
	return [NPOINTS, B-A]

#####SOFA
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





