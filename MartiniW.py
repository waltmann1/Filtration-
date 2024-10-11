from MonomerAbs import MonomerAbs
import numpy as np


class MartiniW(MonomerAbs):

    def __init__(self, n=1, ab=2):

        super(MartiniW, self).__init__()

        self.n = n
        self.ab = ab
        self.length = 2 + self.n + 4 * self.ab
        self.mass = 100 + self.n * 34 + self.ab * 187 + 31


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
        p.bond_names.append("peg")
        p.position.append(np.add(p.position[p.num_beads], [0, 0, factor * .33]))
        p.angle_names.append("mma_join")
        p.angles.append([spot, p.num_beads, p.num_beads+1])
        p.charge.append(0)
        p.mass.append(34)
        p.type.append('EO')
        p.body.append(-1)

        for i in range(1, self.n):
            p.bonds.append([p.num_beads + i, p.num_beads + 1 + i])
            p.bond_names.append("peg")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
            p.angle_names.append("peg")
            p.angles.append([p.num_beads + i - 1, p.num_beads + i, p.num_beads+i + 1])
            p.charge.append(0)
            p.mass.append(34)
            p.type.append('EO')
            p.body.append(-1)

        p.bonds.append([p.num_beads + self.n, p.num_beads + 1 + self.n])
        p.bond_names.append("peg")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
        p.angle_names.append("peg")
        p.angles.append([self.n + p.num_beads - 1, self.n + p.num_beads, self.n + p.num_beads + 1])
        p.charge.append(0)
        p.mass.append(58)
        p.type.append('Na')
        p.body.append(-1)

        p.bonds.append([p.num_beads + 1 + self.n, p.num_beads + 2 + self.n])
        p.bond_names.append("SC1-Na")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
        p.angle_names.append("mma_join")
        p.angles.append([self.n + p.num_beads, self.n + p.num_beads + 1, self.n + p.num_beads + 2])
        p.charge.append(0)
        p.mass.append(28)
        p.type.append('SC1')
        p.body.append(-1)

        p.bonds.append([p.num_beads + 2 + self.n, p.num_beads + 3 + self.n])
        p.bond_names.append("SC1-Na")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
        p.angle_names.append("Na-flex-Na")
        p.angles.append([self.n + p.num_beads + 1, self.n + p.num_beads + 2, self.n + p.num_beads + 3])
        p.charge.append(0)
        p.mass.append(58)
        p.type.append('Na')
        p.body.append(-1)

        p.bonds.append([p.num_beads + 3 + self.n, p.num_beads + 4 + self.n])
        p.bond_names.append("SC1-Na")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
        p.angle_names.append("mma_join")
        p.angles.append([self.n + p.num_beads + 2, self.n + p.num_beads + 3, self.n + p.num_beads + 4])
        p.charge.append(1)
        p.mass.append(43)
        p.type.append('SQd')
        p.body.append(-1)

        for i in range(1, self.ab):
            p.bonds.append([p.num_beads + self.n + i * 4, p.num_beads + 1 + self.n + i * 4])
            p.bond_names.append("SC1-Na")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
            p.angle_names.append("peg")
            p.angles.append([self.n + p.num_beads - 1 + i * 4, self.n + p.num_beads + i * 4,
                             self.n + p.num_beads + 1 + i * 4])
            p.charge.append(0)
            p.mass.append(58)
            p.type.append('Na')
            p.body.append(-1)

            p.bonds.append([p.num_beads + 1 + self.n + i * 4, p.num_beads + 2 + self.n + i * 4])
            p.bond_names.append("SC1-Na")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
            p.angle_names.append("mma_join")
            p.angles.append([self.n + p.num_beads + i * 4, self.n + p.num_beads + 1 + i * 4,
                             self.n + p.num_beads + 2 + i * 4])
            p.charge.append(0)
            p.mass.append(28)
            p.type.append('SC1')
            p.body.append(-1)

            p.bonds.append([p.num_beads + 2 + self.n + i * 4, p.num_beads + 3 + self.n + i * 4])
            p.bond_names.append("SC1-Na")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
            p.angle_names.append("Na-flex-Na")
            p.angles.append([self.n + p.num_beads + 1 + i * 4, self.n + p.num_beads + 2 + i * 4,
                             self.n + p.num_beads + 3 + i * 4])
            p.charge.append(0)
            p.mass.append(58)
            p.type.append('Na')
            p.body.append(-1)

            p.bonds.append([p.num_beads + 3 + self.n + i * 4, p.num_beads + 4 + self.n + i * 4])
            p.bond_names.append("SC1-Na")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .28]))
            p.angle_names.append("mma_join")
            p.angles.append([self.n + p.num_beads + 2 + i * 4, self.n + p.num_beads + 3 + i * 4,
                             self.n + p.num_beads + 4 + i * 4])
            p.charge.append(1)
            p.mass.append(43)
            p.type.append('SQd')
            p.body.append(-1)

        p.bonds.append([p.num_beads + self.n + self.ab * 4, p.num_beads + 1 + self.n + self.ab * 4])
        p.bond_names.append("SC1-Na")
        p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
        p.angle_names.append("peg")
        p.angles.append([self.n + p.num_beads - 1 + self.ab * 4, self.n + p.num_beads + self.ab * 4,
                         self.n + p.num_beads + 1 + self.ab * 4])
        p.charge.append(0)
        p.mass.append(31)
        p.type.append('SP2')
        p.body.append(-1)
        #p.type[-1] = "SP2" #more polar bead for oh group at the end

        for i in range(p.num_beads, p.num_beads + self.length):
            p.monomer_indexes[spot].append(i)
        p.num_beads += self.length