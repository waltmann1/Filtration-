from MonomerAbs import MonomerAbs
import numpy as np


class QVP(MonomerAbs):

    def __init__(self, with_ion=True):

        super(QVP, self).__init__()
        self.length = 1
        self.mass = 1
        self.with_ion = with_ion
        if self.with_ion:
            self.length += 1
            self.mass += 1

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

        p.charge.append(1)
        p.mass.append(1)
        p.type.append('QP')
        p.body.append(-1)

        p.monomer_indexes[spot].append(p.num_beads)

        if self.with_ion:
            p.position.append(np.add(p.position[-1], [0, 0, factor * 1.0]))
            p.charge.append(-1)
            p.mass.append(1)
            p.type.append('QMi')
            p.body.append(-1)
            p.monomer_indexes[spot].append(p.num_beads + 1)

        p.num_beads += self.length
