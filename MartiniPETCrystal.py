from __future__ import division
import numpy as np
from AtomisticPET import AtomisticPET
from SequentialPolymer import SequentialPolymer
from numpy import linalg as la
import copy as cp


class MartiniPETCrystal(object):

    def __init__(self, x, y, z, topology_gone=False):

        self.chains = []
        self.direction = 2


        l_v = AtomisticPET().get_lattice_vectors()
        lengths = [la.norm(l_v[i]) for i in range(3)]
        x_num = int(np.floor(x / lengths[0]))
        y_num = int(np.floor(y / lengths[1]))
        z_num = int(np.floor(z / lengths[2]))
        nums = [x_num, y_num, z_num]

        sequence = ['MartiniPET' for _ in range(z_num)]
        poly = SequentialPolymer(sequence, spiral=False, pet=True)
        poly.restrain_completely(topology_gone=topology_gone)

        for i in range(2 * x_num):
            for j in range(2 * y_num):
                current = cp.deepcopy(poly)
                current.align([0,0,1])
                current.shift(np.add(np.multiply(i, l_v[0]), np.multiply(j, l_v[1])))
                if current.position[0][0] < x and current.position[0][1] < y and current.position[0][0] > 0 and current.position[0][1] > 0:
                    self.chains.append(current)


    def get_chains(self):
        return self.chains

    def shift(self, vector):

        for chain in self.chains:
            chain.shift(vector)

    def orient(self, direction):

        for chain in self.chains:
            chain.position = np.array(chain.position)
            temp = cp.deepcopy(chain.position[:, direction])
            chain.position[:, direction] = chain.position[:, self.direction]
            chain.position[:, self.direction] = temp
        self.direction = direction

    def small_Na(self):

        for chain in self.chains:
            for ind, tipe in enumerate(chain.type):
                if tipe == "Na":
                    chain.type[ind] = "SNa"

    def index_by_type(self, type):

        index_name = "MARTINIPETCrystal_" + type
        for chain in self.chains:
            for i in range(len(chain.position)):
                if chain.type[i] == type:
                    if not index_name in chain.indexed[i]:
                        chain.indexed[i].append(index_name)

    def remove_restraints(self):
        for chain in self.chains:
            chain.restraints = []
            chain.restraint_names = []

    def add_restraints(self):

        for chain in self.chains:
            chain.restrain_completely()

    def restrain_underneath(self):

        max_z = 0
        for chain in self.chains:
            #print(max_z, chain.position[0])
            if chain.position[0][2] > max_z:
                max_z = chain.position[0][2]

        for chain in self.chains:
            #print("looking")
            if chain.position[0][2] != max_z:
                chain.restrain_completely()
            else:
                chain.index_all("MARTINIPETCrystal_" + "top_unrestrained")
                #print("restraining")

    def restrain_ends(self):

        for chain in self.chains:
            chain.restrain_ends()




