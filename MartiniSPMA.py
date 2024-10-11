from MonomerAbs import MonomerAbs
import numpy as np


class MartiniSPMA(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniSPMA, self).__init__()
        self.length = 3
        self.mass = 208

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

        p.position.append(np.add(p.position[-1], [0, 0, factor * .5]))
        p.charge.append(0)
        p.mass.append(28)
        p.type.append('SC1')
        p.body.append(-1)

        p.bonds.append([p.num_beads, p.num_beads + 1])
        p.bond_names.append("peg")

        p.angles.append([spot, p.num_beads, p.num_beads + 1])
        p.angle_names.append('mma_join')

        p.position.append(np.add(p.position[-1], [0, 0, factor * .5]))
        p.charge.append(-1)
        p.mass.append(80)
        p.type.append('Qa')
        p.body.append(-1)

        p.bonds.append([p.num_beads + 1, p.num_beads + 2])
        p.bond_names.append("EHMA2")

        p.angles.append([p.num_beads, p.num_beads + 1, p.num_beads + 2])
        p.angle_names.append('peg')

        for i in range(self.length):
            p.monomer_indexes[spot].append(p.num_beads + 1)
        p.num_beads += self.length
