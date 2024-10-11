from MonomerAbs import MonomerAbs
import numpy as np


class PEGMEMA(MonomerAbs):

    def __init__(self, n=6, with_ion=False):

        super(PEGMEMA, self).__init__()
        self.length = n
        self.mass = n

    def add_to_polymer(self, p, spot, up=True):

        factor = 1
        if not up:
            factor = -1

        p.bonds.append([spot, p.num_beads])
        p.bond_names.append('sidechain')
        p.position.append(np.add(p.position[spot], [0, 0, factor * 1.0]))

        for i in range(self.length):
            if i != 0:
                p.bonds.append([p.num_beads -1 + i, p.num_beads + i])
                p.bond_names.append('sidechain')
                p.position.append(np.add(p.position[-1], [0, 0, factor * 1.0]))

            p.charge.append(0)
            p.mass.append(1)
            p.type.append('L')
            p.body.append(-1)

        for i in range(self.length):
            p.monomer_indexes[spot].append(p.num_beads + i)
        p.num_beads += self.length
