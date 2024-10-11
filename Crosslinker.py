from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class Crosslinker(PolyAbs):

    def __init__(self, n, spiral=False, tight=False, peg=4, k=False, pet=False, mma=False, with_ion=False, ps=False, flanker=True):

        self.spiral = spiral
        self.peg=peg
        self.k = k
        self.flanker = flanker
        sequence = ["A"] + ["B"] * n + ["A"]
        super(Crosslinker, self).__init__(sequence, seperation=4 - 3 * tight, pet=pet, mma=mma, with_ion=with_ion, ps=ps)

        self.sequence = sequence
        self.n =n


    def build_chain(self, sequence, mma=False, with_ion=True, ps=False):

        mon_length = .38

        points = self.linear_points(len(sequence), mon_length)
        if self.spiral:
            points = self.spiral_points(len(sequence), arc=mon_length, separation=self.seperation)

        for i in range(len(sequence)):
            self.position.append(points[i])
            self.type.append(sequence[i])
            self.mass.append(50)
            self.charge.append(0)
            self.body.append(-1)
            self.monomer_indexes.append([i])
            if i != 0:
                k_orig = 150000
                k = k_orig / 2.479
                d = 140
                beta = np.sqrt(k / (2 * d))
                self.bond_names.append("go_morse_0.38_" + str(d) + "_" + str(beta))
                self.bonds.append([i-1, i])
            if i !=0 and i !=1:
                k_angle = 10
                theta =180
                self.angle_names.append("go_harm_" + str(theta) + "_" + str(k_angle))
                self.angles.append([i, i-1, i-2])
            self.length += 1
            self.num_beads += 1
        if self.flanker==True:
            start = self.position[0]
            end = self.position[-1]
            end_index = len(self.position) - 1
            self.body[0] = 0
            self.body[-1] = 1
            offsets = [[0,1,1], [0,1,-1], [0, -1, 1], [0, -1, -1]]
            offsets = np.multiply(offsets, .47)
            for i in range(4):
                self.type.append("FL")
                self.position.append(np.add(start, offsets[i]))
                self.mass.append(1)
                self.charge.append(0)
                self.body.append(0)
                self.monomer_indexes[0].append(self.length)
                self.bond_names.append("go_harm_0.47_1500")
                self.bonds.append([self.length, 0])
                self.length += 1
                self.num_beads += 1



            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length -1, 0, self.length -3])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 1, 0, self.length - 2])

            self.angle_names.append("go_harm_180_10000")
            self.angles.append([self.length - 1, 0, self.length - 4])

            self.angle_names.append("go_harm_180_10000")
            self.angles.append([self.length - 3, 0, self.length - 2])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 3, 0, self.length - 4])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 2, 0, self.length - 4])

            for i in range(4):
                self.type.append("FL")
                self.position.append(np.add(end, offsets[i]))
                self.mass.append(1)
                self.charge.append(0)
                self.body.append(1)
                self.monomer_indexes[end_index].append(self.length)
                self.bond_names.append("go_harm_0.47_1500")
                self.bonds.append([self.length, end_index])
                self.length += 1
                self.num_beads += 1

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 1, end_index, self.length - 3])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 1, end_index, self.length - 2])

            self.angle_names.append("go_harm_180_10000")
            self.angles.append([self.length - 1, end_index, self.length - 4])

            self.angle_names.append("go_harm_180_10000")
            self.angles.append([self.length - 3, end_index, self.length - 2])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 3, end_index, self.length - 4])

            self.angle_names.append("go_harm_90_10000")
            self.angles.append([self.length - 2, end_index, self.length - 4])




