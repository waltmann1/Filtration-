from MonomerAbs import MonomerAbs
import numpy as np


class MartiniPEOPEMAN(MonomerAbs):

    def __init__(self, N, with_ion=False):

        super(MartiniPEOPEMAN, self).__init__(with_ion=False)

        self.n = N
        self.length = 4 + self.n
        self.mass = 100 + 34 * self.n + 72 + 5


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


        p.bonds.append([p.num_beads, p.num_beads+1])
        p.bond_names.append("peo")
        p.position.append(np.add(p.position[p.num_beads], [0, 0, factor * .33]))
        p.angle_names.append("mma_join")
        p.angles.append([spot, p.num_beads, p.num_beads+1])
        p.charge.append(0)
        p.mass.append(34)
        p.type.append('EO')
        p.body.append(-1)


        for i in range(1, self.n):
            p.bonds.append([p.num_beads + i, p.num_beads + 1 + i])
            p.bond_names.append("peo")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
            p.angle_names.append("peo")
            p.angles.append([p.num_beads + i - 1, p.num_beads + i, p.num_beads+i + 1])
            p.charge.append(0)
            p.mass.append(34)
            p.type.append('EO')
            p.body.append(-1)


        p.constraints.append([p.num_beads + self.n, p.num_beads + self.n + 1])
        p.constraint_names.append("SCY-STY")


        p.position.append(np.add(p.position[p.num_beads + self.n], [0, 0, factor * .22]))
        p.charge.append(0)
        p.mass.append(26)
        p.angle_names.append('peo')
        p.angles.append([p.num_beads, p.num_beads + self.n, p.num_beads + self.n + 1])
        p.type.append('STY')
        p.body.append(-1)
        #p.monomer_indexes[spot].append(p.num_beads + self.n + 1)


        p.constraints.append([p.num_beads + self.n + 1, p.num_beads + 2 + self.n])
        p.constraint_names.append('STY-STY')

        p.position.append(np.add(p.position[p.num_beads + 1 + self.n], [0, -.27 * np.sin(np.deg2rad(30)), factor * .27 * np.cos(np.deg2rad(30))]))
        p.charge.append(0)
        p.mass.append(26)
        p.type.append('STY')
        p.body.append(-1)
        #p.monomer_indexes[spot].append(p.num_beads + 3)
        p.angles.append([p.num_beads + self.n, p.num_beads + 1 + self.n, p.num_beads + 2 + self.n])
        p.angle_names.append("B-R-R")

        p.constraints.append([p.num_beads + 2 + self.n, p.num_beads + 3 + self.n])
        p.constraint_names.append('STY-STY')

        p.position.append(
            np.add(p.position[p.num_beads + 1 + self.n], [0, .27 * np.sin(np.deg2rad(30)), factor * .27 * np.cos(np.deg2rad(30))]))
        p.charge.append(0)
        p.mass.append(26)
        p.type.append('STY')
        p.body.append(-1)
        #p.monomer_indexes[spot].append(p.num_beads + 4)

        p.constraints.append([p.num_beads + 1 + self.n, p.num_beads + 3 + self.n])
        p.constraint_names.append('STY-STY')

        for i in range(p.num_beads, p.num_beads + self.length):
            p.monomer_indexes[spot].append(i)
        p.num_beads += self.length
