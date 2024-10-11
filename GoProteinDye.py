from DyeAbs import DyeAbs
import numpy as np
from MartiniProtein import MartiniProtein


class GoProteinDye(DyeAbs):

    def __init__(self, PDBname, itp):


        self.mp = MartiniProtein(PDBname, itp, top=True)
        super(GoProteinDye, self).__init__(index=0)

    def build(self, with_ion=False):

        for i in range(len(self.mp.backbone_atoms)):
            self.position.append(self.mp.position[self.mp.backbone_indices[i] - 1])
            if self.mp.backbone_atoms[i][3] == "LYS" and i !=278:
                self.type.append("A")
            else:
                self.type.append("BB")
            self.mass.append(self.mp.mass[self.mp.backbone_indices[i] - 1])
            self.charge.append(0)

        self.position = np.array(self.position)
        self.position = np.subtract(self.position, np.average(self.position, axis=0))
        self.position = list(self.position)


    def turn_off_reactivity(self):

        for ind, t in enumerate(self.type):
            if t == "A":
                self.type[ind] = "BB"

    def add_topology(self):

        self.index = -1


        for bond in self.mp.backbone_bonds:
            r = bond[3]
            k = bond[4]

            k = k / 2.479

            d = 123

            beta = np.sqrt(k / (2 * d))

            self.bond_names.append("go_morse_" +str(r) + "_" + str(d) + "_" + str(beta))

            #if  self.mp.backbone_indices.index(bond[1]) < 100:
            self.bonds.append([self.mp.backbone_indices.index(bond[0]), self.mp.backbone_indices.index(bond[1])])
        #"""
        for ind, bond in enumerate(self.mp.backbone_disulfide_bonds):
            r = bond[3]
            k = bond[4]
            k = k/2.479
            d = 91

            beta = np.sqrt(k / (2 * d))
            #print(r, d, beta)
            self.bond_names.append("go_morse_" + str(r) + "_" +  str(d) + "_" + str(beta))
            self.bonds.append([self.mp.backbone_indices.index(bond[0]), self.mp.backbone_indices.index(bond[1])])

        for ind, bond in enumerate(self.mp.contact_pairs):
            c6 = bond[3]
            c12 = bond[4]


            sigma = (c6/c12) ** (-1/6)

            eps = c12 / (sigma**12) / 4
            eps = eps / 2.479

            self.bond_names.append("go_lj_" + str(sigma) + "_" + str(eps))
            self.bonds.append([self.mp.backbone_indices.index(bond[0]), self.mp.backbone_indices.index(bond[1])])
        #"""

        for angle in self.mp.backbone_angles:
            #print(angle)

            one = self.mp.backbone_indices.index(angle[0])
            two = self.mp.backbone_indices.index(angle[1])
            three = self.mp.backbone_indices.index(angle[2])

            theta = angle[4]
            k = angle[5]

            k = k/2.479

            self.angle_names.append("go_harm_" + str(theta) + "_" + str(k))
            self.angles.append([one, two, three])


        #quit()
