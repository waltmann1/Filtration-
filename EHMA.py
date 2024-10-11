from MonomerAbs import MonomerAbs
import numpy as np

class EHMA(MonomerAbs):

    def __init__(self, with_ion=False):

        super(EHMA, self).__init__()
        self.length = 2
        self.mass = 2

    def add_to_polymer(self, p, spot, up=True):

        p.bonds.append([spot, p.num_beads])
        p.bond_names.append('sidechain')

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[spot], [0, 0, factor * 1.0]))

        p.charge.append(0)
        p.mass.append(1)
        p.type.append('B')
        p.body.append(-1)

        p.bonds.append([p.num_beads, p.num_beads + 1])
        p.bond_names.append('sidechain')

        p.position.append(np.add(p.position[-1], [0, 0, factor * 1.0]))

        p.charge.append(0)
        p.mass.append(1)
        p.type.append('B')
        p.body.append(-1)

        p.monomer_indexes[spot].append(p.num_beads)
        p.monomer_indexes[spot].append(p.num_beads + 1)
        p.num_beads += self.length
