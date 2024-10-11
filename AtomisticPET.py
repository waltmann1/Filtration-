
from MonomerAbs import MonomerAbs
import numpy as np
import numpy.linalg as la

class AtomisticPET(MonomerAbs):

    def __init__(self):

        super(AtomisticPET, self).__init__()
        self.length = 22
        self.mass = 184 * 2

    def add_to_polymer(self, p, spot, up=True):

        positions, h_positions = self.get_pos()
        start = p.num_beads
        if start:
            start -= 1
            positions = np.add(positions, positions[6])
            h_positions = np.add(h_positions, positions[6])
            positions = np.add(positions, p.positions[start])
            h_positions = np.add(h_positions, p.positions[start])
            p.bonds.append(start, start+1)
            p.bond_names.append("C-C")
        else:
            start -= 1

        p.positions.append(positions[6])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)

        p.positions.append(h_positions[3])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 1, start + 2)
        p.bond_names.append("C-H")

        p.positions.append(h_positions[2])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start+1, start + 3)
        p.bond_names.append("C-H")

        p.positions.append(positions[5])
        p.type.append("O")
        p.mass.append(16)
        p.charge.append(0)
        p.bonds.append(start+1, start + 4)
        p.bond_names.append("C-O")

        p.positions.append(positions[1])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 5, start + 4)
        p.bond_names.append("C-O")

        p.positions.append(positions[4])
        p.type.append("O")
        p.mass.append(16)
        p.charge.append(0)
        p.bonds.append(start + 6, start + 5)
        p.bond_names.append("C-Odouble")

        p.positions.append(positions[0])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 7, start + 5)
        p.bond_names.append("C-C")

        p.positions.append(positions[2])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 7, start + 8)
        p.bond_names.append("C-Cbenz")

        p.positions.append(h_positions[0])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 8, start + 9)
        p.bond_names.append("C-H")

        p.positions.append(positions[3])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 7, start + 10)
        p.bond_names.append("C-Cbenz")

        p.positions.append(h_positions[1])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 10, start + 11)
        p.bond_names.append("C-H")

        p.positions.append(positions[9])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 10, start + 12)
        p.bond_names.append("C-Cbenz")

        p.positions.append(h_positions[4])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 12, start + 13)
        p.bond_names.append("C-H")

        p.positions.append(positions[10])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 8, start + 14)
        p.bond_names.append("C-Cbenz")

        p.positions.append(h_positions[5])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 14, start + 15)
        p.bond_names.append("C-H")

        p.positions.append(positions[7])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 14, start + 16)
        p.bond_names.append("C-Cbenz")
        p.bonds.append(start + 12, start + 16)
        p.bond_names.append("C-Cbenz")

        p.positions.append(positions[8])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 16, start + 17)
        p.bond_names.append("C-C")

        p.positions.append(positions[12])
        p.type.append("O")
        p.mass.append(16)
        p.charge.append(0)
        p.bonds.append(start + 17, start + 18)
        p.bond_names.append("C-Odouble")

        p.positions.append(positions[11])
        p.type.append("O")
        p.mass.append(16)
        p.charge.append(0)
        p.bonds.append(start + 17, start + 19)
        p.bond_names.append("C-O")

        p.positions.append(positions[13])
        p.type.append("C")
        p.mass.append(12)
        p.charge.append(0)
        p.bonds.append(start + 19, start + 22)
        p.bond_names.append("C-O")

        p.positions.append(h_positions[6])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 22, start + 20)
        p.bond_names.append("C-H")

        p.positions.append(h_positions[7])
        p.type.append("H")
        p.mass.append(1)
        p.charge.append(0)
        p.bonds.append(start + 20, start + 21)
        p.bond_names.append("C-H")
        p.num_beads += 22

        for i in range(start, start + self.length):
            p.monomer_indexes[spot].append(p.num_beads)

    def get_pos(self):

        opos = [[0, -0.078, 0.371], [0, -0.181, 0.237], [0.100, 0.170, 0.438], [-0.117, -0.257, 0.433],
                [-0.098, -0.415, 0.170], [0.052, 0.000, 0.175],
                [0.045, -0.085, 0.040]]

        ho_pos = [[0.180, 0.312, 0.385], [-0.204, -0.468, 0.382], [0.290, 0.111, 0.111], [-0.246, -0.137, 0.000]]

        lpos = [self.zero_zero_half_inversion(pos) for pos in opos]
        #print("lpos", lpos)
        hlpos = [self.zero_zero_half_inversion(pos) for pos in opos]

        opos.extend(lpos)
        opos = [self.fractional_to_atomic(pos) for pos in opos]
        ho_pos.extend(hlpos)
        ho_pos = [self.fractional_to_atomic(pos) for pos in hlpos]

        return opos, ho_pos


    def zero_zero_half_inversion(self, vector):
        return [-vector[0], -vector[1], 1 - vector[2]]

    def fractional_to_atomic(self, vector):
        #lattice parmaeters in  nanometers
        a = .456
        b = .594
        c = 1.075
        xy = -0.404
        yz = -0.356
        xz = -0.429
        x_lattice_vector = [1, 0, 0]
        x_lattice_vector = np.multiply(a / la.norm(x_lattice_vector), x_lattice_vector)
        y_lattice_vector = [xy, 1, 0]
        y_lattice_vector = np.multiply(b / la.norm(y_lattice_vector), y_lattice_vector)
        z_lattice_vector = [xz, yz, 1]
        z_lattice_vector = np.multiply(c / la.norm(z_lattice_vector), z_lattice_vector)

        lattice_vectors = [x_lattice_vector, y_lattice_vector, z_lattice_vector]

        atomic = np.sum([np.multiply(lattice_vectors[i], vector[i]) for i in range(3)], axis=0)
        return atomic

    def get_lattice_vectors(self, ):
        a = .456
        b = .594
        c = 1.075
        xy = -0.404
        yz = -0.356
        xz = -0.429
        x_lattice_vector = [1, 0, 0]
        x_lattice_vector = np.multiply(a / la.norm(x_lattice_vector), x_lattice_vector)
        y_lattice_vector = [xy, 1, 0]
        y_lattice_vector = np.multiply(b / la.norm(y_lattice_vector), y_lattice_vector)
        z_lattice_vector = [xz, yz, 1]
        z_lattice_vector = np.multiply(c / la.norm(z_lattice_vector), z_lattice_vector)

        lattice_vectors = [x_lattice_vector, y_lattice_vector, z_lattice_vector]

        return lattice_vectors


    """
C(1) -0.013 (0.000) -0.088 (-0.078) 0.369 (0.371)
C(2) -0.028 (0.000) -0.181 (-0.181) 0.230 (0.237)
C(3) 0.101 (0.100) 0.176 (0,170) 0.436 (0.438)
C(4) -0.115 (-0.117) -0.263 (-0.257) 0.434 (0.433)
H(8) 0.180 0.312 0.385
H(9) -0.204 -0.468 0.382
Non-group atoms
0(6) -0.099 (-0.098) -0.4t9 (-0.415) 0.169 (0:170)
0(7) 0.054 (0.052) 0.001 (0.000) 0.173 (0.175)
C(5) 0.045 (0.045) -0.090 (-0.085) 0.038 (0.040)
H(10)~ ,0.290 0.111 0,111
H(11):~ -0.246 -0.137 0.000 
https://link.springer.com/content/pdf/10.1007/BF01202554.pdf
a = 4.56
b = 5.94
c = 10.75
alpha = 98.5
beta = 112
gamma = 118
xy = -0.404
yz = -0.356
xy = -0.429
    """
