import numpy as np
from numpy import linalg as la


class DyeAbs(object):

        def __init__(self, index=0):
            self.center = [0, 0, 0]
            self.mass = []
            self.position = []
            self.charge = []
            self.type = []
            self.orientation = [1, 0, 0, 0]
            self.index = index
            self.bonds = []
            self.bond_names = []
            self.angles = []
            self.angle_names = []
            self.build()

        def build(self):

            print("I am an abstract class and I ain't building shit")

        def shift(self, vector):

            for ind, site in enumerate(self.position):
                self.position[ind] = np.add(site, vector)

            self.center = np.add(self.center, vector)

        def align(self, quat):

            for ind, site in enumerate(self.position):
                self.position[ind] = quat.orient(site)

            self.orientation = quat.q

            return [quat.q[3], quat.q[0], quat.q[1], quat.q[2]]

        def max_radius(self, counterions=False):

            #print("calculating max radius")
            cen = self.rigid_center_of_mass()
            #print("center", cen)
            max = 0
            for ind, pos in enumerate(self.position):
                dist = la.norm(np.subtract(pos, cen))
                #print(cen, pos, dist)
                if dist > max and self.type[ind][-1] != 'i':
                    max = dist

            return max

        def calculate_inertia_tensor(self):

            tensor = np.zeros((3, 3))
            position_array = np.subtract(self.position, self.center)
            mass_array = self.mass
            for idx in range(len(mass_array)):
                tensor[0][0] += (np.square(position_array[idx][1]) + np.square(position_array[idx][2])) * mass_array[
                    idx]
                tensor[1][1] += (np.square(position_array[idx][2]) + np.square(position_array[idx][0])) * mass_array[
                    idx]

                tensor[2][2] += (np.square(position_array[idx][0]) + np.square(position_array[idx][1])) * mass_array[
                    idx]
                tensor[0][1] -= position_array[idx][0] * position_array[idx][1] * mass_array[idx]
                tensor[0][2] -= position_array[idx][0] * position_array[idx][2] * mass_array[idx]
                tensor[1][2] -= position_array[idx][1] * position_array[idx][2] * mass_array[idx]
                tensor[1][0] = tensor[0][1]
                tensor[2][0] = tensor[0][2]
                tensor[2][1] = tensor[1][2]

                values, vectors = la.eig(tensor)
            return values

        def rigid_center_of_mass(self):

            length = len(self.type)
            index = len(self.type) - 1
            while self.type[index][-1] == "i":
                length -= 1
                index -= 1
            mass = self.mass[:length]
            mass_array = np.array(mass)
            position_array = np.array(self.position[:length])
            array = [0,0,0]

            for i in range(length):
                for x in range(3):
                    array[x] += mass_array[i] * position_array[i][x]

            return array / np.sum(mass_array)

        def enforce_cubic_bc(self, box_length):

            self.image = [[0, 0, 0] for _ in range(len(self.position))]
            half = box_length / 2
            changed = False
            for ind1, position in enumerate(self.position):
                if not self.is_ion(ind1):
                    for ind2 in range(0, 3):
                        if position[ind2] > box_length + half:
                            raise ValueError("dye bead is greater than a box length outside")
                        elif position[ind2] > half:
                            self.position[ind1][ind2] -= box_length
                            self.image[ind1][ind2] += 1
                            changed = True
                            #print("fixed", self.position[ind1], self.image[ind1])
                        elif position[ind2] < -1 * (box_length + half):
                            raise ValueError("dye bead is greater than a box length outside")
                        elif position[ind2] < -1 * half:
                            self.position[ind1][ind2] += box_length
                            self.image[ind1][ind2] -= 1
                            changed = True
            return changed

        def is_ion(self, ind):

            return self.type[ind][-1] == "i"

