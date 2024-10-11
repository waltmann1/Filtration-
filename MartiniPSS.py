from MonomerAbs import MonomerAbs
import numpy as np
from MartiniPS import MartiniPS


class MartiniPSS(MartiniPS):

    def __init__(self, with_ion=False):

        super(MartiniPSS, self).__init__()
        self.mass += 72

    def add_to_polymer(self, p, spot, up=True):

        super(MartiniPSS, self).add_to_polymer(p, spot, up=up)

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[spot], [0, 0, factor * (.22 + 2 * .27 * np.cos(np.deg2rad(30)))]))
        p.charge.append(-1)
        p.mass.append(72)
        p.type.append('Qa')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)

        p.bonds.append([p.num_beads, p.num_beads -1])
        p.bond_names.append("STY-STY")
        p.bonds.append([p.num_beads, p.num_beads - 2])
        p.bond_names.append("STY-STY")
        p.dihedrals.append([p.num_beads, p.num_beads - 1, p.num_beads - 2, p.num_beads - 3])
        p.dihedral_names.append("Felipe")
        p.num_beads += 1
        self.length += 1