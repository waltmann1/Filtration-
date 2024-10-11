from __future__ import division
import numpy as np
from PolyAbs import PolyAbs

class SequentialPolymer(PolyAbs):

    def __init__(self, sequence, spiral=True, tight=False, peg=4, k=False, pet=False, mma=False, with_ion=True, ps=False):

        self.spiral = spiral
        self.peg=peg
        self.k = k
        super(SequentialPolymer, self).__init__(sequence, seperation=4 - 3 * tight, pet=pet, mma=mma, with_ion=with_ion, ps=ps)
        self.sequence = sequence

    def build_chain(self, sequence, mma=False, with_ion=True, ps=False):

        mon_length = .84
        if mma:
            mon_length = .289
        elif ps:
            mon_length = .25


        points = self.linear_points(len(sequence), mon_length)
        if self.spiral:
            points = self.spiral_points(len(sequence), arc=mon_length, separation=self.seperation)

        for i in range(len(sequence)):
            self.position.append(points[i])
            self.type.append('B')
            self.mass.append(1)
            self.charge.append(0)
            self.body.append(-1)
            self.monomer_indexes.append([i])
            if i != 0:
                self.bond_names.append('polybond')
                if mma:
                    self.bond_names[-1] = "MMABackbone"
                elif ps:
                    self.bond_names[-1] = "SCY-SCY"
                self.bonds.append([i-1, i])
            if self.k or mma or ps:
                if i !=0 and i !=1:
                    self.angle_names.append('k'+ str(self.k))
                    if mma:
                        self.angle_names[-1] = "MMABackbone"
                    elif ps:
                        self.angle_names[-1] = "B-B-B"
                    self.angles.append([i-2, i-1, i])
            self.length += 1
            self.num_beads += 1
        for spot, thing in enumerate(sequence):
            imp = __import__(thing, fromlist=[''])
            if thing == "PEGMEMA":
                mon = getattr(imp, thing)(n=self.peg, with_ion=with_ion)
            else:
                print(imp, thing)
                mon = getattr(imp, thing)(with_ion=with_ion)
            mon.add_to_polymer(self, spot, up=(spot % 2 == 0))
        if ps:
            for mon_index, mon in enumerate(self.monomer_indexes[:-1]):
                self.angles.append([mon[1], mon[0], self.monomer_indexes[mon_index + 1][0]])
                self.angle_names.append("R-B-B")

    def martini_build(self, sequence):

        for i in range(len(sequence)):
            self.monomer_indexes.append([])

        for spot, thing in enumerate(sequence):
            imp = __import__(thing, fromlist=[''])
            mon = getattr(imp, thing)()
            mon.add_to_polymer(self, spot, up=(spot % 2 == 0))
