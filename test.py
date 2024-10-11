from  __future__ import division
import sys
sys.path.append('/home/cwj8781/Filtration')

from SequentialPolymer import SequentialPolymer
from OneBeadDye import OneBeadDye
from Solution import Solution
from MartiniProtein import MartiniProtein
from SequenceGenerator import SequenceGenerator
#from Simulation import Simulation
import copy as cp
import numpy as np
from numpy import linalg as la

box_length =50

polys = []
#protein = MartiniProtein("protein_CG.pdb", "PETase") 
#protein.shift(-protein.rigid_center_of_mass())
#protein.shift([box_length/2,box_length/2,box_length/2])
#print(protein.itp)
#polys.append(protein)

mono_names = ['MartiniMMA', 'MartiniPEGMA']
mono_probabilities = [0.5, 0.5]
num_chains = 4
chain_length = 100

rando_generator = SequenceGenerator(num_chains, chain_length, mono_names, mono_probabilities)



for ind, sequence in enumerate(rando_generator.get_sequences()):
	poly = SequentialPolymer(sequence, spiral=False, mma=True)
	alignment = (np.random.random_sample((3,)) - .5)
	alignment = np.divide(alignment, la.norm(alignment))
	poly.align([0,0,1])
	poly.shift(-1 * poly.geometric_center())
	position = [ind*10, 15, 0 ] 
	poly.shift(position)
	polys.append(poly)

for poly in polys:
	poly.enforce_cubic_bc(box_length)
	poly.shift([box_length/2, box_length/2, box_length/2])
	



d = OneBeadDye()
d.shift([0, 5, 0])

dyes = [d]

s = Solution(dyes, polys, box_length)

s.martini_init("pets")

s.dump_gsd()
