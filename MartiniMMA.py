from MonomerAbs import MonomerAbs
import numpy as np


class MartiniMMA(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniMMA, self).__init__()
        self.length = 1
        self.mass = 100

    def add_to_polymer(self, p, spot, up=True):

        p.type[spot] = "SC1"
        p.mass[spot] = 41

        p.bonds.append([spot, p.num_beads])
        p.bond_names.append("SC1-Na")

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[spot], [0, 0, factor * .28]))
        p.charge.append(0)
        p.mass.append(59)
        p.type.append('Na')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)
        p.num_beads += self.length