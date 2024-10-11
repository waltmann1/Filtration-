from __future__ import division
import numpy as np
from numpy import linalg as la
from Quaternion import QuaternionBetween
import math



class PolyAbs(object):

    def __init__(self, sequence, index=-1, seperation=4, pet=False, mma=False, with_ion=True, ps=False):


        self.position = []
        self.type = []
        self.bonds = []
        self.bond_names = []
        self.angles = []
        self.angle_names = []
        self.charge = []
        self.length = 0
        self.rigid_count = 0
        self.mass = []
        self.index = index
        self.sequence = sequence
        self.monomer_indexes = []
        self.num_beads = 0
        self.body = []
        self.seperation = seperation
        self.dihedrals = []
        self.dihedral_names = []
        self.constraints = []
        self.constraint_names = []
        self.restraints = []
        self.restraint_names = []


        if pet:
            self.martini_build(sequence)
        else:
            self.build_chain(sequence, mma=mma, with_ion=with_ion, ps=ps)
        self.image = []
        self.itp = False
        self.indexed = [[] for _ in range(len(self.position))]

    def build_chain(self, sequence, mma=False, with_ion=True):
        print("I'm an abstract class and I ain't building shit")

    def martini_build(self, sequence):
        print("I literally don't know how to do a martini build, I am PolyAbs")

    def align(self, vec):
        q = QuaternionBetween(self.chain_vector(), vec)
        for x in range(len(self.position)):
            self.position[x] = q.orient(self.position[x])

    def align_to_q(self, q):
        for x in range(len(self.position)):
            self.position[x] = q.orient(self.position[x])

    def shift(self, vector):

        for ind, site in enumerate(self.position):
            self.position[ind] = np.add(site, vector)

    def chain_vector(self):

        return np.subtract(self.position[-1], self.position[0])

    def rigid_center_of_mass(self):

        mass = self.mass

        mass_array = np.array(mass)
        position_array = np.array(self.position)
        print(mass_array.shape)
        print(position_array.shape)
        return np.sum(self.pos_x_mass(position_array, mass_array), axis=0) / np.sum(mass_array)

    def pos_x_mass(self, position_array, mass_array):

        y = np.zeros_like(position_array)
        for ind in range(len(mass_array)):
            y[ind][0] = position_array[ind][0] * mass_array[ind]
            y[ind][1] = position_array[ind][1] * mass_array[ind]
            y[ind][2] = position_array[ind][2] * mass_array[ind]
        return y

    def moment_inertia(self):
        mass = self.mass

        mass_array = np.array(mass)
        position_array = np.array(self.position)

        cen = self.center_of_mass_arrays(position_array, mass_array)
        position_array = np.subtract(position_array, cen)
        return self.sum_over_xyz(self.pos_x_mass(self.pos_squared(position_array), mass_array))

    def pos_squared(self, position_array):

        y = np.zeros_like(position_array)

        for ind in range(len(position_array)):
            y[ind][0] = position_array[ind][0] * position_array[ind][0]
            y[ind][1] = position_array[ind][1] * position_array[ind][1]
            y[ind][2] = position_array[ind][2] * position_array[ind][2]
        return y

    def sum_over_xyz(self, array):

        final = np.array([0, 0, 0])

        for list in array:
            final[0] += list[0]
            final[1] += list[1]
            final[2] += list[2]
        return final

    def center_of_mass_arrays(self, position_array, mass_array):
        return np.sum(self.pos_x_mass(position_array, mass_array), axis=0) / np.sum(mass_array)

    def max_radius(self):

        cen = self.rigid_center_of_mass()
        max = 0
        for pos in self.position:
            dist = la.norm(np.subtract(pos, cen))
            if dist > max:
                max = dist
        return max

    def spiral_points(self,n, arc=.5, separation=4):
        """generate points on an Archimedes' spiral
        with `arc` giving the length of arc between two points
        and `separation` giving the distance between consecutive
        turnings
        - approximate arc length with circle arc at given distance
        - use a spiral equation r = b * phi
        """

        def p2c(r, phi):
            """polar to cartesian
            """
            return [r * math.cos(phi), r * math.sin(phi), 0]

        # yield a point at origin
        points=  [[0,0,0]]

        # initialize the next point in the required distance
        r = arc
        b = separation / (2 * math.pi)
        # find the first phi to satisfy distance of `arc` to the second point
        phi = float(r) / b
        count = 0
        while count < n - 1:
            points.append(p2c(r, phi))
            # advance the variables
            # calculate phi that will give desired arc length at current radius
            # (approximating with circle)
            phi += float(arc) / r
            r = b * phi
            count += 1
        return points

    def linear_points(self, number, spacing):

        return [[x * spacing, 0, 0] for x in range(number)]

    def geometric_center(self):

        return np.average(self.position, axis=0)

    def enforce_cubic_bc(self, box_length):

        self.image = [[0, 0, 0] for _ in range(len(self.position))]
        half = box_length / 2
        for ind1, position in enumerate(self.position):
            if not self.is_ion(ind1):
                for ind2 in range(0, 3):
                    if position[ind2] > box_length + half:
                        raise ValueError("Polymer bead is greater than a box length outside")
                    elif position[ind2] > half:
                        self.position[ind1][ind2] -= box_length
                        self.image[ind1][ind2] += 1
                    elif position[ind2] < -1 * (box_length + half):
                        raise ValueError("Polymer bead is greater than a box length outside")
                    elif position[ind2] < -1 * half:
                        self.position[ind1][ind2] += box_length
                        self.image[ind1][ind2] -= 1

    def is_ion(self, ind):

        return self.type[ind][-1] == "i"

    def restrain_completely(self, name="normal", topology_gone=False):

        self.restraint_names = []
        self.restraints = []

        for i in range(len(self.position)):
            self.restraints.append(i)
            if topology_gone:
                self.restraint_names.append("times10")
            else:
                self.restraint_names.append(name)
        if topology_gone:
            self.angles = []
            self.angle_names = []
            self.dihedrals = []
            self.dihedral_names = []
            self.bonds = []
            self.bond_names = []
            self.constraint_names = []
            self.constraints = []

    def index_all(self, name):

        for i in range(len(self.position)):
            if name not in self.indexed[i]:
                self.indexed[i].append(name)

    def index_by_type(self, type):

        for i in range(len(self.type)):
            if self.type[i] == type:
                self.indexed[i].append(type)

    def restrain_ends(self, name="normal"):

        for index in self.monomer_indexes[-1]:
            if index not in self.restraints:
                self.restraints.append(index)
                self.restraint_names.append(name)

        for index in self.monomer_indexes[0]:
            if index not in self.restraints:
                self.restraints.append(index)
                self.restraint_names.append(name)

    def itp_bonds_string(self, b_obj):

        string = ""
        for ind, bond in enumerate(self.bonds):
            string += b_obj.martini_string(bond, self.bond_names[ind])
        return string

    def itp_angles_string(self, a_obj):

        string = ""
        for ind, angles in enumerate(self.angles):
            string += a_obj.martini_string(angles, self.angle_names[ind])
        return string

    def itp_constraints_string(self, c_obj):

        string = ""
        for ind, constraint in enumerate(self.constraints):
            string += c_obj.martini_string(constraint, self.constraint_names[ind])
        return string

    def itp_restraints_string(self, p_obj):

        string = ""
        for ind, restraint in enumerate(self.restraints):
            string += p_obj.martini_string(restraint, self.restraint_names[ind])
        return string

    def itp_dihedrals_string(self, d_obj):

        string = ""
        for ind, dihedral in enumerate(self.dihedrals):
            string += d_obj.martini_string(dihedral, self.dihedral_names[ind])
        return string

    def itp_pairs_string(self):

        return "\n"

    def itp_atoms_string(self, name):

        string = ""
        for atom in range(self.num_beads):
            nr = atom + 1
            type = self.type[atom]
            resnum = 1
            residue = name
            atomname = type
            cnr = atom + 1
            charge = self.charge[atom]
            mass = self.mass[atom]
            string += "%s%5s%5d%8s%5s%5d%8.3f%8.3f\n" % (nr, type, resnum, residue, atomname, cnr, charge, mass)
        return string

    def full_itp_string(self, name, a, b, c, d, p):

        str_a = self.itp_angles_string(a)
        str_b = self.itp_bonds_string(b)
        str_c = self.itp_constraints_string(c)
        str_d = self.itp_dihedrals_string(d)
        str_p = self.itp_restraints_string(p)
        str_pairs = self.itp_pairs_string()
        str_atoms = self.itp_atoms_string(name)
        exclusions = np.sum([len(thing) > 0 for thing in [str_a, str_b, str_c, str_d, str_p]])

        string = ""
        string += "[ moleculetype ] \n"
        string += "; name  exclusions\n"
        string += name + "    " + str(exclusions) + "\n"
        string += "\n"
        string += "[ atoms ] \n"
        string += "; nr    type      resnumber     residue      atom    cgnr     charge     mass\n"
        string += str_atoms
        string += "\n"
        if len(str_b) > 0:
            string += "[ bonds ] \n"
            string += "; a     b      funct    param1    param2\n"
            string += str_b
            string += "\n"
        if len(str_c) > 0:
            string += "[ constraints ] \n"
            string += "; a    b    funct   param1\n"
            string += str_c
            string += "\n"
        if len(str_a) > 0:
            string += "[ angles ] \n"
            string += "; a     b      c   param1    param2\n"
            string += str_a
            string += "\n"
        if len(str_d) > 0:
            string += "[ dihedrals ] \n"
            string += "; a     b   c    d  param1    param2\n"
            string += str_d
        if len(str_p) > 0:
            string += "\n[ position_restraints ] \n"
            string += "; a  funct param1    param2 param3\n"
            string += str_p
        if len(str_pairs) > 0:
            string += "\n[ pairs ] \n"
            string += "; a     b      funct    param1    param2\n"
            string += str_pairs

        return string