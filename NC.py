from MonomerAbs import MonomerAbs
import numpy as np


class NC(MonomerAbs):

    def __init__(self, with_ion=True):

        super(NC, self).__init__()
        self.length = 0
        self.mass = 0
        self.with_ion = with_ion
        if self.with_ion:
            self.length = 1
            self.mass = 1

    def add_to_polymer(self, p, spot, up=True):

        p.charge[spot] = -1
        p.type[spot] = 'QM'
        p.mass[spot] = 1

        factor = 1
        if not up:
            factor = -1

        if self.with_ion:
            p.position.append(np.add(p.position[spot], [0, 0, factor * 1.0]))
            p.charge.append(1)
            p.mass.append(1)
            p.type.append('QPi')
            p.body.append(-1)

        if self.with_ion:
            p.monomer_indexes[spot].append(p.num_beads)
        p.num_beads += self.length