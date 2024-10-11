import numpy as np
from GoProteinDye import GoProteinDye
import copy as cp
import random


class FunctionalizedGPD(object):

    def __init__(self, PDBname, itp, chains):


        self.gpd = GoProteinDye(PDBname, itp)
        self.gpd.add_topology()
        self.chains = [cp.deepcopy(chain) for chain in chains]
        self.setup_chains()

    def setup_chains(self):

        link_points = []

        for ind, t in enumerate(self.gpd.type):
            if t == "A":
                link_points.append(ind)


        random.shuffle(link_points)
        self.link_points = link_points[:len(self.chains)]

        center = self.gpd.rigid_center_of_mass()
        offset_magnitude = 10

        for ind, point in enumerate(self.link_points):
            vec = np.subtract(self.gpd.position[point], center)
            self.chains[ind].align(vec)
            offset = np.multiply(np.divide(vec, np.linalg.norm(vec)), offset_magnitude)
            shift = np.add(vec, offset)
            self.chains[ind].shift(shift)

        self.gpd.turn_off_reactivity()

    def shift(self, vec):

        self.gpd.shift(vec)
        for chain in self.chains:
            chain.shift(vec)

    def enforce_cubic_bc(self, box_length):

        self.gpd.enforce_cubic_bc(box_length)
        for chain in self.chains:
            chain.enforce_cubic_bc(box_length)