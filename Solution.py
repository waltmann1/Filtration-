from __future__ import division
import numpy as np
from numpy import linalg as la
import hoomd

import copy as cp
from SequentialPolymer import SequentialPolymer
from PositiveIon import PositiveIon
from NegativeIon import NegativeIon
from DivalentNegativeIon import DivalentNegativeIon
from DivalentPositiveIon import DivalentPositiveIon
import random

class Solution(object):

    def __init__(self, dyes, chains, box_length, functional_bonds=None):

        self.dyes = [cp.deepcopy(dye) for dye in dyes]
        self.chains = [cp.deepcopy(chain) for chain in chains]
        self.dump_context = None
        self.rigid_count = 0
        for dye in self.dyes:
            if dye.index != -1:
                self.rigid_count +=1
        self.num_particles = np.sum([len(dye.mass) for dye in dyes]) + np.sum([len(chain.mass) for chain in self.chains])
        self.dyes = self.reindex()
        self.name = "assembly"
        self.box_length = int(box_length)

        self.o_list = None
        self.index_names = []
        self.index_indices = []
        self.functional_bonds=functional_bonds


    def reindex(self):
        tag = 0
        new_list = []
        added = []
        for index, dye in enumerate(self.dyes):
            if dye.index != -1:
                dye.index = tag
                tag += 1
                new_list.append(dye)
                added.append(index)
        for index, dye in enumerate(self.dyes):
            if index not in added:
                new_list.append(dye)
        return new_list

    def functionalize_dye(self, snap, chain_index, index_in_dye, dye_index, bond_name, all_bond_types):

        bond_number = snap.bonds.N + 1
        snap.bonds.resize(bond_number)
        func_index = np.sum([len(self.dyes[i].position) for i in range(dye_index)]) + index_in_dye + self.rigid_count
        index_for_chain = np.sum([len(self.dyes[i].position) for i in range(len(self.dyes))]) + self.rigid_count + \
                          np.sum([len(self.chains[i].position) for i in range(chain_index)])
        snap.bonds.group[bond_number - 1] = [func_index, index_for_chain]
        snap.bonds.typeid[bond_number - 1] = all_bond_types.index(bond_name)


    def create_system(self):
        """

        :return: system object
        """

        import hoomd
        if self.dump_context is None:
            self.dump_context = hoomd.context.initialize("")
        b_types = self.parse_bond_types()
        p_types = self.parse_particle_names() + ["center"]
        a_types = self.parse_angle_types()
        l = self.box_length

        print(self.num_particles, self.rigid_count, p_types, b_types, a_types, l)
        snap = hoomd.data.make_snapshot(int(self.num_particles + self.rigid_count), particle_types=p_types,
                                        bond_types=b_types, angle_types=a_types,
                                        box=hoomd.data.boxdim(L=l))
        snap.bonds.resize(0)

        for x in range(self.rigid_count):
            snap.particles.position[x] = self.dyes[x].rigid_center_of_mass()
            snap.particles.mass[x] = np.sum([self.dyes[x].mass])
            snap.particles.typeid[x] = p_types.index("center")
            snap.particles.body[x] = x
            snap.particles.moment_inertia[x] = self.dyes[x].calculate_inertia_tensor()

        tag = self.rigid_count
        for dye in self.dyes:
            for x in range(len(dye.bonds)):
                bond_number = snap.bonds.N + 1
                snap.bonds.resize(bond_number)
                snap.bonds.group[bond_number - 1] = np.add(dye.bonds[x], tag)
                snap.bonds.typeid[bond_number - 1] = b_types.index(dye.bond_names[x])

            for x in range(len(dye.angles)):
                angle_number = snap.angles.N + 1
                snap.angles.resize(angle_number)
                snap.angles.group[angle_number - 1] = np.add(dye.angles[x], tag)
                snap.angles.typeid[angle_number - 1] = a_types.index(dye.angle_names[x])
            for x in range(len(dye.position)):
                snap.particles.position[x + tag] = dye.position[x]
                snap.particles.mass[x + tag] = dye.mass[x]
                snap.particles.typeid[x + tag] = p_types.index(dye.type[x])
                if dye.type[x][-1] != 'i':
                    snap.particles.body[x + tag] = dye.index
                else:
                    snap.particles.body[x + tag] = -1
                snap.particles.charge[x + tag] = dye.charge[x]
            tag += len(dye.position)

        for chain in self.chains:
            for x in range(len(chain.bonds)):
                bond_number = snap.bonds.N + 1
                snap.bonds.resize(bond_number)
                snap.bonds.group[bond_number - 1] = np.add(chain.bonds[x], tag)
                snap.bonds.typeid[bond_number - 1] = b_types.index(chain.bond_names[x])

            for x in range(len(chain.angles)):
                angle_number = snap.angles.N + 1
                snap.angles.resize(angle_number)
                snap.angles.group[angle_number - 1] = np.add(chain.angles[x], tag)
                snap.angles.typeid[angle_number - 1] = a_types.index(chain.angle_names[x])

            for x in range(chain.num_beads):
                snap.particles.position[x + tag] = chain.position[x]
                snap.particles.mass[x + tag] = chain.mass[x]
                snap.particles.typeid[x + tag] = p_types.index(chain.type[x])
                snap.particles.body[x + tag] = chain.index
                snap.particles.charge[x + tag] = chain.charge[x]
                print(chain.image[x])
                print(snap.particles.image[x + tag])
                snap.particles.image[x + tag] = chain.image[x]
            tag += chain.num_beads

        if self.functional_bonds is not None:
            for bond in self.functional_bonds:
                chain_index, index_in_dye, dye_index = bond
                bond_name = "go_attach"
                self.functionalize_dye(snap, chain_index, index_in_dye, dye_index, bond_name, b_types)

        sys = hoomd.init.read_snapshot(snap)
        self.dump_map()

        return sys

    def dec3(self, string):

        for ind, char in enumerate(string):
            if char == ".":
                return string[:ind+4]
        return string + ".000"

    def martini_init(self, title):

        self.coordinate_file_chains_only(title)
        self.topology_file_chains_only(title)
        self.dump_xyz(title)

    def dump_xyz(self, title):

        num = np.sum([len(chain.position) for chain in self.chains])

        f = open(title + ".xyz", "w")
        f.write(str(num))
        f.write("\n\n")
        count = 0
        mon_count = 0

        for ind, chain in enumerate(self.chains):
            print(ind, chain.itp)
            for mon in range(chain.num_beads):
                count += 1
                mon_count += 1
                s = "%5s%8.3f%8.3f%8.3f\n" % (
                    chain.type[mon], chain.position[mon][0], chain.position[mon][1], chain.position[mon][2])
                f.write(s)
            mon_count = 0
        f.write(str(self.box_length) + " " + str(self.box_length) + " " + str(self.box_length))
        f.close()

    def coordinate_file_chains_only(self, title):

        num = np.sum([len(chain.position) for chain in self.chains])

        f = open(title + ".gro", "w")
        f.write(title + "\n")
        f.write("%5d\n" % num)
        count = 0
        mon_count = 0

        for ind, chain in enumerate(self.chains):
            resnumber = 0
            for mon in range(chain.num_beads):
                count += 1
                mon_count += 1
                resname = "Po" + str(int(ind+1))
                if chain.type[mon] == "BB" or mon == 0:
                    resnumber += 1
                if chain.itp:
                    resname = chain.sequence[resnumber - 1]
                s = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (
                    resnumber, resname, chain.type[mon], count, chain.position[mon][0], chain.position[mon][1],
                    chain.position[mon][2], 0, 0, 0)
                if len(chain.indexed[mon]) > 0:
                    for index_name in chain.indexed[mon]:
                        #print(index_name)
                        if index_name in self.index_names:
                            self.index_indices[self.index_names.index(index_name)].append(count)
                        else:
                            self.index_names.append(index_name)
                            self.index_indices.append([count])
                f.write(s)
            mon_count = 0
        s = "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n" % \
            (self.box_length, self.box_length, self.box_length, 0, 0, 0, 0, 0, 0)
            #(box[X][X], box[Y][Y], box[Z][Z], box[X][Y], box[X][Z], box[Y][X], box[Y][Z], box[Z][X], box[Z][Y])
        f.write(s)
        f.close()

    def topology_file_chains_only(self, title):

        from Angles import Angles
        from Bonds import Bonds
        from Constraints import Constraints
        from Dihedrals import Dihedrals
        from PositionRestraints import PositionRestraints

        a = Angles()
        b = Bonds()
        c = Constraints()
        d = Dihedrals()
        p = PositionRestraints()

        f = open(title + ".top", 'w')
        itp_names = ["polymers.itp"]
        for chain in self.chains:
            if chain.itp:
                itp_names.append(chain.itp + ".itp")
        itp_string = ""
        f2 = open(itp_names[0], 'w')
        f.write('#include "martini_v2.2P_PS.itp"' + "\n")
        f.write('#include "martini_v2.0_ions.itp"' + "\n")
        f.write('#include "martini_v2.0_solvents.itp"' + "\n")
        for itp_name in itp_names:
            f.write('#include "' + itp_name + '"' + "\n")
        f.write("\n")
        f.write("[ system ]" + "\n")
        f.write(title + "\n")

        f.write("[ molecules ]\n")
        for ind in range(len(self.chains)):
            if not self.chains[ind].itp:
                chain_name = "Po" + str(ind + 1)
                itp_string += "\n"
                itp_string += self.chains[ind].full_itp_string(chain_name, a, b, c, d, p)
            else:
                chain_name = self.chains[ind].itp
            f.write(chain_name + " " + "1" + "\n")
        f2.write(itp_string)
        f.close()
        f2.close()

    def parse_bond_types(self):

        a_types = ["go_attach"]
        for chain in self.chains:
            for c in chain.bond_names:
                if c not in a_types:
                    a_types.append(c)
        for chain in self.dyes:
            for c in chain.bond_names:
                if c not in a_types:
                    a_types.append(c)
        return a_types

    def parse_angle_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.angle_names:
                if c not in a_types:
                    a_types.append(c)
        for chain in self.dyes:
            for c in chain.angle_names:
                if c not in a_types:
                    a_types.append(c)

        return a_types

    def parse_particle_names(self):
        a_types = []
        for chain in self.chains:
            for c in chain.type:
                if c not in a_types:
                    a_types.append(c)
        for dye in self.dyes:
            for d in dye.type:
                if d not in a_types:
                    a_types.append(d)

        return a_types

    def geometric_center(self):

        pos = [np.mean(c.position, axis=0) for c in self.chains]
        pos_dyes = [dye.center for dye in self.dyes]
        pos = pos + pos_dyes
        weights = [len(c.position) for c in self.chains]
        weights_mers = [len(dye.position) for dye in self.dyes]
        weights = weights + weights_mers
        total = np.sum(weights)
        weights = [float(weight)/float(total) for weight in weights]
        center = np.array([0, 0, 0])
        for ind, p in enumerate(pos):
            #print(weights[ind])
            center = np.add(center, np.multiply(p, weights[ind]))
        return center

    def max_radius(self):

        t = [-1 * c for c in self.geometric_center()]
        self.shift(t)
        rs = [la.norm(pos) for chain in self.chains for pos in chain.position]
        self.shift([-1 * g for g in t])
        return np.max(rs)

    def shift(self, vector):

        for dye in self.dyes:
            dye.shift(vector)

        for chain in self.chains:
            chain.shift(vector)

    def center_of_mass(self, center_index):
        for dye in self.dyes:
            if dye.index == center_index:
                return dye.rigid_center_of_mass()
        for chain in self.chains:
            if chain.index == center_index:
                return chain.rigid_center_of_mass()

    def total_rigid_mass(self, center_index):
        for dye in self.dyes:
            if dye.index == center_index:
                return dye.total_rigid_mass()
        for chain in self.chains:
            if chain.index == center_index:
                return np.sum(chain.mass)

    def moment_inertia(self, center_index):
        for dye in self.dye:
            if dye.index == center_index:
                return dye.calculate_interia_tensor()

    def dump_gsd(self, dump_file=None):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """

        filename = dump_file
        if dump_file is None:
            filename = self.name + '.gsd'
        elif dump_file[-4:] != '.gsd':
            filename = dump_file + '.gsd'

        sys = self.create_system()
        #res_map = self.create_res_map(sys)

        hoomd.dump.gsd(filename=filename, period=None, group=hoomd.group.all(), overwrite=True)
        return filename

    def orient_quaternion(self, q):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))
        for chain in self.chains:
            for x in range(chain.num_particles):
                chain.position[x] = q.orient(chain.position[x])
        self.shift(temp)

    def center_at_origin(self):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))

    def dump_map(self, dump_file=None):

        qpi_list = []
        qmi_list = []
        filename = dump_file
        if dump_file is None:
            filename = self.name + '.map'
        elif dump_file[-4:] != '.map':
            filename = dump_file + '.map'
        f_out = open(filename, 'w')
        tag = self.rigid_count
        for dye in self.dyes:
            string = "dye "
            for i in range(len(dye.position)):
                string += str(tag)
                string += " "
                tag += 1
            string += "\n"
            f_out.write(string)
        for chain in self.chains:
            if "QP" in chain.type:
                string = "charged_polymer "
            else:
                string = "rando_polymer "
            for ind, name in enumerate(chain.sequence):
                string += name
                string += " "
                for index in chain.monomer_indexes[ind]:
                    try:
                        chain.type[index][-1] != "i"
                    except IndexError:
                        print(index, chain.type, len(chain.type))
                    if chain.type[index][-1] != "i":
                        string += str(index + tag)
                        string += " "
                    else:
                        if chain.type[index][1] == "P":
                            qpi_list.append(index+tag)
                        else:
                            qmi_list.append(index+tag)
            string += "\n"
            f_out.write(string)
            tag += chain.num_beads
        string = "qmi "
        for index in qmi_list:
            string += str(index)
            string += " "
        string += "\n"
        f_out.write(string)
        string = "qpi "
        for index in qpi_list:
            string += str(index)
            string += " "
        string += "\n"
        f_out.write(string)

        f_out.close()

    def add_dye(self, polys, dyes, dye_index, shift_vector, box_length):
        reset_dye = cp.deepcopy(dyes[dye_index])
        dyes[dye_index].shift(shift_vector)
        if dyes[dye_index].enforce_cubic_bc(box_length):
            dyes[dye_index] = reset_dye
            return False

        obstacles = polys + [dye for ind, dye in enumerate(dyes) if dye.index <= 0 and
                             ind != dye_index]
        positions = np.array([pos for obstacle in obstacles for pos in obstacle.position])
        centers = np.array([dyes[dye_index].rigid_center_of_mass()
                            for _ in range(positions.shape[0])])
        subs = np.subtract(positions, centers)
        dists = np.sum(np.abs(subs) ** 2, axis=-1) ** (1. / 2)
        cut = np.multiply(np.ones_like(dists), dyes[dye_index].max_radius()+1)
        val = np.sum(np.greater(cut, dists))
        print(val)
        if val == 0:
            return True
        else:
            dyes[dye_index] = reset_dye
            return False

    def martini_volume_fraction(self, frame=0):

        types = ['STY', 'C1']

        diameters = [.27, .43, .47]

        found_diameters = []

        for chain in self.chains:
            for tip in chain.type:
                if tip in types:
                    found_diameters.append(diameters[types.index(tip)])
                else:
                    found_diameters.append(diameters[-1])

        vol = 0
        for d in found_diameters:
            vol += (d/2)**3 * 4/3 * np.pi
        frac = vol / (self.box_length**3)
        return frac

    def dump_index_file(self):

        f = open("index.ndx", "w")

        for i in range(len(self.index_names)):
            print(self.index_names, self.index_names[i])
            f.write("[ " + self.index_names[i] + " ]\n")
            strind = ""
            for ndx in self.index_indices[i]:
                strind += str(ndx) + " "
                if len(strind) > 4080:
                    strind += "\n"
                    f.write(strind)
                    strind = ""
            strind += "\n"
            f.write(strind)


class RandomSolution(Solution):

    def __init__(self, rando_generator, charged_generator, dyes, box_length=50, scatter_ions=True, tight=False, peg=6, k=False,
                 pos_k=False, remove_ions=False, divalent=False, total_compensation=False):

        if k:
            k =int(k)
        if pos_k:
            pos_k = int(pos_k)
        polys = []

        with_ion = not remove_ions

        print("random sol with ion", with_ion)

        for ind, sequence in enumerate(rando_generator.get_sequences() + charged_generator.get_sequences()):
            if ind < rando_generator.num_chains or not pos_k:
                poly = SequentialPolymer(sequence, tight=tight, peg=peg, k=k, with_ion=with_ion)
            else:
                poly = SequentialPolymer(sequence, tight=tight, peg=peg, k=pos_k, with_ion=with_ion)
            alignment = (np.random.random_sample((3,)) - .5)
            alignment = np.divide(alignment, la.norm(alignment))
            poly.align(alignment)
            poly.shift(-1 * poly.geometric_center())
            position = (np.random.random_sample((3,)) - .5) * box_length
            poly.shift(position)


            if scatter_ions:
                for ind, position in enumerate(poly.position):
                    if poly.type[ind][-1] == 'i':
                        poly.position[ind] = (np.random.random_sample((3,)) - .5) * box_length

            poly.enforce_cubic_bc(box_length)
            polys.append(poly)

        if remove_ions:
            poly_charge = np.sum([np.sum(chain.charge) for chain in polys])
            total_negative_charge = -1 * np.sum([np.sum(chain.charge) for chain in polys if np.sum(chain.charge) < 0])
            total_positive_charge = np.sum([np.sum(chain.charge) for chain in polys if np.sum(chain.charge) > 0])
            if poly_charge > 0 and not divalent:
                for _ in range(int(poly_charge)):
                    dyes.append(NegativeIon())
            elif poly_charge > 0 and divalent and not total_compensation:
                for _ in range(int(poly_charge/2)):
                    dyes.append(DivalentNegativeIon())
            elif poly_charge < 0 and divalent and not total_compensation:
                for _ in range(int(poly_charge/2)):
                    dyes.append(DivalentPositiveIon())
            elif divalent and total_compensation:
                for _ in range(int(total_positive_charge/2)):
                    dyes.append(DivalentNegativeIon())
                for _ in range(int(total_negative_charge/2)):
                    dyes.append(DivalentPositiveIon())
                #print(poly_charge, total_positive_charge, total_negative_charge, len(dyes))
            else:
                for _ in range(int(-1 * poly_charge)):
                    dyes.append(PositiveIon())

        #for dye in dyes:
         #   val = 1
         #   while val != 0:
         #      point = (np.random.random_sample((3,)) - .5) * box_length
         #      tests = [dye.max_radius() > la.norm(np.subtract(point, point2)) for poly in polys + dyes for point2 in poly.position]
         #      val = sum(tests)
         #      print(dye.max_radius(), point, val)
         #   dye.shift(point)

        for dye_index, dye in enumerate(dyes):
            shift_vector = (np.random.random_sample((3,)) - .5) * box_length
            dye.shift(shift_vector)
            #while not self.add_dye(polys, dyes, dye_index, shift_vector, box_length):
            #    shift_vector = (np.random.random_sample((3,)) - .5) * box_length
            if scatter_ions:
                for ind, position in enumerate(dye.position):
                    if dye.type[ind][-1] == 'i':
                        dye.position[ind] = (np.random.random_sample((3,)) - .5) * box_length


            dye.enforce_cubic_bc(box_length)

        super(RandomSolution, self).__init__(dyes, polys, box_length)


class LatticeSolution(Solution):

    def __init__(self, rando_generator, charged_generator, dyes, spacing, n, box_length=50, scatter_ions=True, peg=6, k=False,
                 pos_k=False):

        sequences = rando_generator.get_sequences() + charged_generator.get_sequences()
        angle_k = [int(k) for _ in range(rando_generator.num_chains)] + [pos_k for _ in range(charged_generator.num_chains)]
        save_order = cp.deepcopy(sequences)
        np.random.shuffle(sequences)
        angle_k_order = np.zeros(len(sequences))
        for i in range(len(sequences)):
            angle_k_order[i] = int(angle_k[save_order.index(sequences[i])])
        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        z = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy, zz = np.meshgrid(x, y, z)
        pos = []
        polys = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    pos.append([xx[i][j][k], yy[i][j][k], zz[i][j][k]])
        for ind, seq in enumerate(sequences):
            poly = SequentialPolymer(seq, tight=True, peg=peg, k=int(angle_k_order[ind]))
            poly.shift(pos[ind])
            poly.enforce_cubic_bc(box_length)


            if scatter_ions:
                for ind, position in enumerate(poly.position):
                    if poly.type[ind][-1] == 'i':
                        poly.position[ind] = (np.random.random_sample((3,)) - .5) * box_length

            polys.append(cp.deepcopy(poly))

        for dye in dyes:
            dye.shift((np.random.random_sample((3,)) - .5) * box_length)
            if scatter_ions:
                for ind, position in enumerate(dye.position):
                    if dye.type[ind][-1] == 'i':
                        dye.position[ind] = (np.random.random_sample((3,)) - .5) * box_length
            dye.enforce_cubic_bc(box_length)

        super(LatticeSolution, self).__init__(dyes, polys, box_length)


class OrderedSolution(Solution):

    def __init__(self, polys, dyes, spacing, n, box_length=50):


        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        z = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy, zz = np.meshgrid(x, y, z)
        pos = []
        total_list = polys+dyes

        random.shuffle(total_list)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    pos.append([xx[i][j][k], yy[i][j][k], zz[i][j][k]])
        print(len(pos))
        for ind,thing in enumerate(total_list):
            if i < len(polys):
                thing.shift(pos[ind])
                thing.enforce_cubic_bc(box_length)

        super(OrderedSolution, self).__init__(dyes=dyes, chains=polys, box_length=box_length)


class FunctionalizedSolution(Solution):

    def __init__(self, polys, f_dyes, spacing, n, box_length=50):


        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        z = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy, zz = np.meshgrid(x, y, z)
        pos = []

        n_f = int(np.floor(n / np.power(len(f_dyes), .3333)))
        n_f = int(np.ceil(30 / spacing))
        pos_n = []
        total_list = f_dyes + polys

        for i in range(0, n, n_f):
            for j in range(0, n, n_f):
                for k in range(0, n, n_f):
                    pos_n.append([xx[i][j][k], yy[i][j][k], zz[i][j][k]])


        #random.shuffle(total_list)
        index_list = list(range(n))
        for index in index_list:
            if index % n_f == 0:
                index_list.remove(index)

        for i in index_list:
            for j in index_list:
                for k in index_list:
                    temp_pos = [xx[i][j][k], yy[i][j][k], zz[i][j][k]]
                    pos.append(temp_pos)

        random.shuffle(pos)
        random.shuffle(pos_n)
        pos_n.extend(pos)
        pos = pos_n
        print(len(pos))
        for ind,thing in enumerate(total_list):
            #if i < len(polys):
            thing.shift(pos[ind])
            thing.enforce_cubic_bc(box_length)


        new_polys = polys
        new_dyes = []
        functional_bonds= []

        for ind, fdye in enumerate(f_dyes):
            new_dyes.append(cp.deepcopy(fdye.gpd))
            for ind2, link_index in enumerate(fdye.link_points):
                index_in_dye = link_index
                dye_index = ind
                chain_index = len(new_polys)
                functional_bonds.append([chain_index, index_in_dye, dye_index])
                fdye.chains[ind2].enforce_cubic_bc(box_length)
                new_polys.append(fdye.chains[ind2])

        super(FunctionalizedSolution, self).__init__(dyes=new_dyes, chains=new_polys, box_length=box_length,
                                                     functional_bonds=functional_bonds)



class GeneralRandomSolution(Solution):

    def __init__(self, polys, dyes, box_length=50):

        for poly in polys:
            alignment = (np.random.random_sample((3,)) - .5)
            alignment = np.divide(alignment, la.norm(alignment))
            poly.align(alignment)
            poly.shift(-1 * poly.geometric_center())
            position = (np.random.random_sample((3,)) - .5) * box_length
            poly.shift(position)
            poly.enforce_cubic_bc(box_length)

        to_remove = []
        for dye_index, dye in enumerate(dyes):
            shift_vector = (np.random.random_sample((3,)) - .5) * (box_length - 2 * dye.max_radius())
            dye.shift(shift_vector)

            dye.enforce_cubic_bc(box_length)

        super(GeneralRandomSolution, self).__init__(dyes, polys, box_length)


class TwoComponentSolution(GeneralRandomSolution):

    def __init__(self, poly, num_poly, dye, num_dye, box_length=50):

        polys = [cp.deepcopy(poly) for _ in range(num_poly)]
        dyes = [cp.deepcopy(dye) for _ in range(num_dye)]

        super(TwoComponentSolution, self).__init__(polys, dyes, box_length)
