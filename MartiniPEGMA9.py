from MonomerAbs import MonomerAbs
import numpy as np


class MartiniPEGMA9(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniPEGMA9, self).__init__(with_ion=False)

        self.n=9
        self.length = self.n + 1
        self.mass = 100 + self.n * 34


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
        p.type.append('SN0')
        p.body.append(-1)

        for i in range(1, self.n):
            p.bonds.append([p.num_beads + i, p.num_beads + 1 + i])
            p.bond_names.append("peg")
            p.position.append(np.add(p.position[-1], [0, 0, factor * .33]))
            p.angle_names.append("peg")
            p.angles.append([p.num_beads + i - 1, p.num_beads + i, p.num_beads+i + 1])
            p.charge.append(0)
            p.mass.append(34)
            p.type.append('SN0')
            p.body.append(-1)

        #p.type[-1] = "SP2" #more polar bead for oh group at the end

        for i in range(p.num_beads, p.num_beads + self.length):
            p.monomer_indexes[spot].append(i)
        p.num_beads += self.length