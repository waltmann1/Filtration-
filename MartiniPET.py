from __future__ import division
from MonomerAbs import MonomerAbs
import numpy as np
from AtomisticPET import AtomisticPET
from numpy import linalg as la

class MartiniPET(MonomerAbs):

    def __init__(self, with_ion=False):

        super(MartiniPET, self).__init__()
        self.length = 5
        self.mass = 180

    def add_to_polymer(self, p, spot, up=True):

        start = p.num_beads
        com = np.array(self.centers_of_mass())

        starting_pos = [0, 0, 0]

        if start:
            com = np.add(com, com[0])
            starting_pos = p.position[start-1]


        com = np.add(com, starting_pos)
        #print("starting pos")
        #print(starting_pos)
        #print("com")
        #print(com)
        p.position.append(com[0])
        p.charge.append(0)
        p.mass.append(58)
        p.type.append('Na')
        #p.type.append('SNa')
        p.body.append(-1)

        if start:
            p.bonds.append([start, start - 1])
            p.bond_names.append('Na-Na')
            p.angles.append([start - 2, start - 1, start])
            p.angle_names.append("STY-Na-Na")

        p.position.append(com[1])
        p.charge.append(0)
        p.mass.append(24)
        p.type.append('STY')
        p.body.append(-1)

        p.bonds.append([start, start + 1])
        p.bond_names.append('Na-STY')

        if start:
            p.angles.append([start - 1, start, start + 1])
            p.angle_names.append("Na-Na-STY")

        p.position.append(com[2])
        p.charge.append(0)
        p.mass.append(24)
        p.type.append('STY')
        p.body.append(-1)

        p.constraints.append([start + 1, start + 2])
        p.constraint_names.append('STY-STY')
        #p.bonds.append([start + 1, start + 2])
        #p.bond_names.append('STY-STY')
        p.angles.append([start, start + 1, start + 2])
        p.angle_names.append("Na-STY-STY")

        p.position.append(com[3])
        p.charge.append(0)
        p.mass.append(24)
        p.type.append('STY')
        p.body.append(-1)

        #p.bonds.append([start + 1, start + 3])
        #p.bond_names.append('STY-STY')
        p.constraints.append([start + 1, start + 3])
        p.constraint_names.append('STY-STY')

        #p.bonds.append([start + 2, start + 3])
        #p.bond_names.append('STY-STY')
        p.constraints.append([start + 2, start + 3])
        p.constraint_names.append('STY-STY')

        #p.angle_names.append("STY-STY-STY")
        #p.angles.append([start + 1, start + 2, start + 3])
        #p.angle_names.append("STY-STY-STY")
        #p.angles.append([start + 2, start + 3, start + 1])
        #p.angle_names.append("STY-STY-STY")
        #p.angles.append([start + 3, start + 1, start + 2])

        p.dihedrals.append([start + 3, start + 2, start + 1, start])
        p.dihedral_names.append("Felipe")

        p.position.append(com[4])
        p.charge.append(0)
        p.mass.append(58)
        p.type.append('Na')
        #p.type.append('SNa')
        p.body.append(-1)

        p.bonds.append([start + 2, start + 4])
        p.bond_names.append('STY-Na_1')

        p.bonds.append([start + 3, start + 4])
        p.bond_names.append('STY-Na_2')

        p.angles.append([start + 1, start + 2, start + 4])
        p.angle_names.append("STY-STY-Na")

        p.dihedrals.append([start + 4, start + 3, start + 2, start + 1])
        p.dihedral_names.append("Felipe")

        for i in range(start, start + self.length):
            p.monomer_indexes[spot].append(p.num_beads)
        p.num_beads += self.length

    def centers_of_mass(self):

        atomic_pos, h_pos = AtomisticPET().get_pos()
        Na1_indices = [1, 4, 5, 6]
        Na1_pos = [atomic_pos[index] for index in Na1_indices]
        Na1_mass = [12, 16, 16, 14]
        Na1_mass_frac = np.true_divide(Na1_mass, np.sum(Na1_mass))
        Na1_center_of_mass = np.sum([Na1_mass_frac[i] * Na1_pos[i] for i in range(len(Na1_mass))], axis=0)

        Sty1_center_of_mass = atomic_pos[0]
        Sty2_center_of_mass = np.add(np.multiply(atomic_pos[10], 12/13), np.multiply(h_pos[5], 1/13))
        Sty3_center_of_mass = np.add(np.multiply(atomic_pos[9], 12 / 13), np.multiply(h_pos[4], 1 / 13))

        """
        print("0 2",la.norm(np.subtract(atomic_pos[0], atomic_pos[2])))
        print("2 10",la.norm(np.subtract(atomic_pos[2], atomic_pos[10])))
        print("10 7",la.norm(np.subtract(atomic_pos[10], atomic_pos[7])))
        print("7 9",la.norm(np.subtract(atomic_pos[7], atomic_pos[9])))
        print("9 3",la.norm(np.subtract(atomic_pos[3], atomic_pos[9])))
        print("3 0",la.norm(np.subtract(atomic_pos[0], atomic_pos[3])))
        print("3 0", la.norm(np.subtract(atomic_pos[0], atomic_pos[3])))
        """
        #print("sty na", la.norm(np.subtract(Sty1_center_of_mass, Na1_center_of_mass)))


        Na2_indices = [8, 11, 12, 13]
        Na2_pos = [atomic_pos[index] for index in Na2_indices]
        Na2_mass = [12, 16, 16, 14]
        Na2_mass_frac = np.true_divide(Na2_mass, np.sum(Na2_mass))
        Na2_center_of_mass = np.sum([Na2_mass_frac[i] * Na2_pos[i] for i in range(len(Na2_mass))], axis=0)

        #print("na na", 2 * la.norm(np.subtract([0,0,0], Na1_center_of_mass)))
        #print("na sty", la.norm(np.subtract(Na2_center_of_mass, Sty3_center_of_mass)))
        #print("na sty", la.norm(np.subtract(Na2_center_of_mass, Sty2_center_of_mass)))

        return [Na1_center_of_mass, Sty1_center_of_mass, Sty2_center_of_mass, Sty3_center_of_mass, Na2_center_of_mass]

    def angles(self):

        names = ["Na-Na-STY", "Na-STY-STY", "STY-STY-STY", "STY-STY-Na", "STY-Na-Na"]

        ap = AtomisticPET()
        next_vec = ap.fractional_to_atomic([0, 0, 1])
        last_vec = -next_vec

        com = self.centers_of_mass()

        last_na = np.add(last_vec, com[4])
        na_na_vec1 = np.subtract(com[0], last_vec)
        na_sty_vec = np.subtract(com[1], com[0])
        dot = np.sum([na_na_vec1[i] * na_sty_vec[i] for i in range(3)])
        angle0 = np.rad2deg(np.arccos(dot / la.norm(na_sty_vec) / la.norm(na_na_vec1)))


        sty_sty_vec = np.subtract(com[2], com[1])
        dot = np.sum([na_sty_vec[i] * sty_sty_vec[i] for i in range(3)])
        angle1 = np.rad2deg(np.arccos(dot / la.norm(na_sty_vec) / la.norm(sty_sty_vec)))

        angle2 = 120
        sty_na_vec = np.subtract(com[4], com[2])
        dot = np.sum([sty_na_vec[i] * sty_sty_vec[i] for i in range(3)])
        angle3 = np.rad2deg(np.arccos(dot / la.norm(sty_na_vec) / la.norm(sty_sty_vec)))

        next_na = np.add(com[0], next_vec)
        na_na_vec2 = np.subtract(next_na, com[4])
        dot = np.sum([sty_na_vec[i] * na_na_vec2[i] for i in range(3)])
        angle4 = np.rad2deg(np.arccos(dot / la.norm(sty_na_vec) / la.norm(na_na_vec2)))

        angles = [angle0, angle1, angle2, angle3, angle4]

        return angles, names

    def bonds(self):

        names = ["Na-Na", "Na-STY", "STY-Na_1", "STY-Na_2"]

        ap = AtomisticPET()
        next_vec = ap.fractional_to_atomic([0, 0, 1])
        last_vec = -next_vec

        com = self.centers_of_mass()

        last_na = np.add(last_vec, com[4])
        bond_1 = la.norm(np.subtract(last_na, com[0]))
        bond_2 = la.norm(np.subtract(com[0], com[1]))
        bond_3 = la.norm(np.subtract(com[2], com[4]))
        bond_4 = la.norm(np.subtract(com[3], com[4]))

        bonds = [bond_1, bond_2, bond_3, bond_4]

        return bonds, names