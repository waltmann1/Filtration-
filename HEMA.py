from MonomerAbs import MonomerAbs
import numpy as np


class HEMA(MonomerAbs):

    def __init__(self, with_ion=False):

        super(HEMA, self).__init__()
        self.length = 1
        self.mass = 1

    def add_to_polymer(self, p, spot, up=True):
        p.type[spot] = 'B'

        factor = 1
        if not up:
            factor = -1

        p.bonds.append([spot, p.num_beads])
        p.bond_names.append('sidechain')

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[spot], [0, 0, factor * 1.0]))

        p.charge.append(0)
        p.mass.append(1)
        p.type.append('L')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)

        p.num_beads += self.length