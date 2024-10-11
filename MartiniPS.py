from MonomerAbs import MonomerAbs
import numpy as np


class MartiniPS(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniPS, self).__init__()
        self.length = 3
        self.mass = 119

    def add_to_polymer(self, p, spot, up=True):

        p.type[spot] = "SCY"
        p.mass[spot] = 41

        p.constraints.append([spot, p.num_beads])
        p.constraint_names.append("SCY-STY")

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[spot], [0, 0, factor * .22]))
        p.charge.append(0)
        p.mass.append(26)
        p.type.append('STY')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)

        p.constraints.append([p.num_beads, p.num_beads + 1])
        p.constraint_names.append('STY-STY')

        factor = 1
        if not up:
            factor = -1

        p.position.append(np.add(p.position[p.num_beads], [0, -.27 * np.sin(np.deg2rad(30)), factor * .27 * np.cos(np.deg2rad(30))]))
        p.charge.append(0)
        p.mass.append(26)
        p.type.append('STY')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)
        p.angles.append([spot, p.num_beads, p.num_beads + 1])
        p.angle_names.append("B-R-R")


        p.constraints.append([p.num_beads + 1, p.num_beads + 2])
        p.constraint_names.append('STY-STY')


        p.position.append(
            np.add(p.position[p.num_beads], [0, .27 * np.sin(np.deg2rad(30)), factor * .27 * np.cos(np.deg2rad(30))]))
        p.charge.append(0)
        p.mass.append(26)
        p.type.append('STY')
        p.body.append(-1)
        p.monomer_indexes[spot].append(p.num_beads)

        p.constraints.append([p.num_beads, p.num_beads + 2])
        p.constraint_names.append('STY-STY')

        p.num_beads += self.length

