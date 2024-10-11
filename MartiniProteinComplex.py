from __future__ import division
import numpy as np
from MartiniProtein import MartiniProtein
from Quaternion import QuaternionBetween

class MartiniProteinComplex(object):

    #assumes one pdb file but multiple itps
    #

    def __init__(self, PDBs, itps):

        self.polys = []
        if len(PDBs) != len(itps):
            raise ValueError("need same number of pdbs as itps")
        for i in range(len(PDBs)):
            self.polys.append(MartiniProtein(PDBs[i], itps[i]))

    def shift(self, vector):

        for poly in self.polys:
            poly.shift(vector)


    def get_chains(self):

        return self.polys

    def rigid_center_of_mass(self):

        masses = []
        positions = []
        for poly in self.polys:
            masses.append(np.sum(poly.mass))
            positions.append(poly.rigid_center_of_mass())
        mass_array = np.array(masses)
        position_array = np.array(positions)
        #print(mass_array)
        #print(position_array)
        return np.sum(self.polys[0].pos_x_mass(position_array, mass_array), axis=0) / np.sum(mass_array)

    def index_all(self, name):
        for poly in self.polys:
            poly.index_all(name)

    def align(self, vec):

        chain1_vec = self.polys[0].chain_vector()
        q = QuaternionBetween(chain1_vec, vec)

        for poly in self.polys:
            temp = poly.rigid_center_of_mass()
            poly.shift(-temp)
            poly.align_to_q(q)
            poly.shift(temp)
