from MonomerAbs import MonomerAbs
import numpy as np


class MartiniHEMA(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniHEMA, self).__init__()
        self.length = 2
        self.mass = 131

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

        p.bonds.append([p.num_beads, p.num_beads + 1])
        p.bond_names.append("SC1-Na")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
        p.angle_names.append("mma_join")
        p.angles.append([spot, p.num_beads, p.num_beads + 1])
        p.charge.append(0)
        p.mass.append(31)
        p.type.append('SP2')
        p.body.append(-1)


        p.monomer_indexes[spot].append(p.num_beads)
        p.monomer_indexes[spot].append(p.num_beads + 1)
        p.num_beads += self.length