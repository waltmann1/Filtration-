from MonomerAbs import MonomerAbs
import numpy as np


class OldSPMMA(MonomerAbs):

    def __init__(self, with_ion=True):

        super(OldSPMMA, self).__init__()
        self.length = 2
        self.mass = 2
        self.with_ion = with_ion
        if with_ion:
            self.length = 3
            self.mass = 3

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
        p.charge.append(-1)
        p.mass.append(1)
        p.type.append('QM')
        p.body.append(-1)

        if self.with_ion:
            p.position.append(np.add(p.position[-1], [0, 0, factor * 0.8]))
            p.charge.append(1)
            p.mass.append(1)
            p.type.append('QPi')
            p.body.append(-1)

        p.monomer_indexes[spot].append(p.num_beads)
        p.monomer_indexes[spot].append(p.num_beads + 1)
        if self.with_ion:
            p.monomer_indexes[spot].append(p.num_beads + 2)
        p.num_beads += self.length