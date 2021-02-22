import sys
import numpy as np
import libsofaTerminal as sofa
import scipy.sparse as sp
import time

#####SOFA
sofa.open("Armadillo.scn")

ROOT=sofa.root()

def setRes(R):
	GRID=ROOT.getObject("SparseGrid")
	GRID.clear()
	GRID.getData("cellWidth").setValue(R)
	GRID.init()

	MAPPING=ROOT.getObject("DeformArmadillo1/mapping")
	MAPPING.init()
	
	CONTAINER=ROOT.getObject("DeformArmadillo1/Container")
	POS=CONTAINER.getData("position").getValue(False)
	
	STATE=ROOT.getObject("DeformArmadillo1/mstate")
	STATE.resize(0)
	STATE.resize(CONTAINER.getData("nbPoints").getValue())
	STATE.getData("rest_position").setValue(POS)
	STATE.getData("reset_position").setValue(POS)
	STATE.reset()

	MASS=ROOT.getObject("DeformArmadillo1/mass")
	MASS.getData("indices").setValue(CONTAINER.getData("points").getValue())
	
	ROOT.getContext("DeformArmadillo1").propagateTopologicalChanges()
	
	#MASS.handleTopologyChange	
	#FEM=ROOT.getObject("DeformArmadillo1/FEM")
	#FEM.handleTopologyChange()

sofa.gui()

A=time.perf_counter()
for i in range(1, 100):
	sofa.step()
B=time.perf_counter()
print("TIME =", B-A)

setRes(1.2)
sofa.gui()

A=time.perf_counter()
for i in range(1, 100):
	sofa.step()
B=time.perf_counter()
print("TIME =", B-A)



