from __future__ import division
import numpy as np
from PolyAbs import PolyAbs
import copy as cp

class MartiniProtein(PolyAbs):

    def __init__(self, PDBname, itp, top=False):

        #the PDB should be from martinize.py
        sequence = PDBname
        super(MartiniProtein, self).__init__(sequence, seperation=4 - 3 * False, pet=False, mma=False, with_ion=False)
        self.itp = itp
        self.petase_active_site = [117, 118, 266, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315,
                                            316, 360, 361, 423, 424, 425, 426, 427, 428]
        self.petase_active_site_backbone = [117, 266, 271, 312, 360, 423, 427]

        self.top = top
        if top:
            self.get_topology(itp)

    def build_chain(self, sequence, mma=False, with_ion=False, ps=False):

        pdb_name = sequence
        self.sequence = []

        f = open(pdb_name, 'r')

        data = f.readlines()
        res_number = -.77
        total_count = 0
        for line in data:
            s = line.split()
            new_res = False
            if s[0] == "ATOM":
                #print(s)
                type = s[2]
                res_name = s[3]
                res_column = 4
                if s[res_column].isalpha():
                    res_column =5
                if int(s[res_column]) != res_number:
                    new_res = True
                res_number = int(s[res_column])
                x = float(s[res_column + 1]) / 10
                y = float(s[res_column + 2]) / 10
                z = float(s[res_column + 3]) / 10
                self.position.append([x, y, z])
                self.type.append(type)
                self.mass.append(40)
                self.charge.append(0)
                self.body.append(-1)
                if new_res:
                    self.monomer_indexes.append([total_count])
                    self.sequence.append(res_name)
                else:
                    self.monomer_indexes[-1].append(total_count)
                total_count += 1
                self.num_beads += 1


    def get_topology(self, itp_name):

        f = open(itp_name, "r")

        data = f.readlines()

        backbone_atoms = []
        backbone_indices = []

        self.backbone_bonds = []
        self.other_bonds = []

        self.backbone_angles = []
        self.other_angles = []

        self.backbone_dihedrals = []
        self.other_dihedrals = []

        self.disulfide_bonds = []
        self.backbone_disulfide_bonds = []

        self.contact_pairs = []

        atoms =False
        bonds = False
        bb_bonds = False
        angles= False
        bb_angles = False
        dihedrals = False
        bb_dihedrals=False
        constraints=False
        pairs = False

        for line in data:
            s = line.split()
            if len(s) > 1:
                if s[0] == "[":
                    atoms = False
                    bonds = False
                    bb_bonds = False
                    angles = False
                    bb_angles = False
                    dihedrals = False
                    bb_dihedrals = False
                    constraints = False
                    pairs=False
                if s[1] == "atoms":
                    atoms=True
                elif s[1] == "bonds":
                    bonds=True
                elif s[1] == "angles":
                    angles=True
                elif s[1] == "dihedrals":
                    dihedrals=True
                elif s[1] == "pairs":
                    pairs=True
                elif s[1] == "constraints":
                    constraints = True
                elif s[1] == "Backbone" and bonds:
                    bb_bonds =True
                elif s[1] == "Backbone" and angles:
                    bb_angles =True
                elif s[1] == "Backbone" and dihedrals:
                    bb_dihedrals =True
                elif s[1] == "Sidechain" and bonds:
                    bb_bonds =False
                elif s[1] == "Sidechain" and angles:
                    bb_angles =False
                elif s[1] == "Sidechain" and dihedrals:
                    bb_dihedrals =False

            if len(s) > 4:
                if atoms and s[4] == "BB":
                    backbone_atoms.append([int(s[0]), str(s[1]), int(s[2]), str(s[3]), str(s[4]), int(s[5]), float(s[6])])
                    backbone_indices.append(int(s[0]))
                elif bonds and bb_bonds:
                    self.backbone_bonds.append([int(s[0]), int(s[1]), int(s[2]), float(s[3]), float(s[4])])
                elif constraints and len(s) == 7:
                    self.disulfide_bonds.append([int(s[0]), int(s[1]), int(s[2]), float(s[3])])
                elif bonds and not bb_bonds:
                    self.other_bonds.append([int(s[0]), int(s[1]), int(s[2]), float(s[3]), float(s[4])])
                elif angles and bb_angles:
                    self.backbone_angles.append([int(s[0]), int(s[1]), int(s[2]), int(s[3]), float(s[4]), float(s[5])])
                elif angles and not bb_angles:
                    self.other_angles.append([int(s[0]), int(s[1]), int(s[2]), int(s[3]), float(s[4]), float(s[5])])
                elif dihedrals and bb_dihedrals:
                    self.backbone_dihedrals.append([int(s[0]), int(s[1]), int(s[2]), int(s[3]), int(s[4]), float(s[5]), float(s[6])])
                elif dihedrals and not bb_dihedrals:
                    self.other_dihedrals.append([int(s[0]), int(s[1]), int(s[2]), int(s[3]), int(s[4]), float(s[5]), float(s[6])])
                elif pairs:
                    self.contact_pairs.append([int(s[0]), int(s[1]), int(s[2]), float(s[3]), float(s[4])])

        self.backbone_atoms = backbone_atoms
        self.backbone_indices = backbone_indices
        self.backbone_disulfide_bonds = self.get_backbone_disulfide()

    def get_backbone_disulfide(self):

        new_pairs = []
        for di in self.disulfide_bonds:
            first = di[0]
            while first not in self.backbone_indices:
                first -= 1
            second = di[1]
            while second not in self.backbone_indices:
                second -= 1
            first_pos = self.position[first - 1]
            second_pos = self.position[second - 1]
            dist = np.linalg.norm(np.subtract(first_pos, second_pos))
            new_pairs.append([first,second, 1, dist, 90000])
        return  new_pairs


    def martini_build(self, sequence):

        print("Never do this")



    def index_backbone_cut_terminals(self, start, end):

        self.index_by_type("BB")
        the_list = []
        for spot, indexed in enumerate(self.indexed):
            if "BB" in indexed:
                the_list.append(spot)
                self.indexed[spot] = indexed[:-1]
        if end != 0:
            the_list = the_list[start:-end]
        else:
            the_list = the_list[start:]

        for thing in the_list:
            self.indexed[thing].append("BB_cut_" + str(start) + "_" + str(end))
            print(start, end, self.indexed[thing])

    def index_petase_active_site_backbone(self):

        for index in self.petase_active_site_backbone:
            self.indexed[index].append("PETase_Active_Site_Backbone")

    def index_petase_active_site(self):

        for index in self.petase_active_site:
            self.indexed[index].append("PETase_Active_Site")

    def go_coordinate_file(self, title, length=None):

        num = len(self.backbone_atoms)

        f = open(title + ".gro", "w")
        f.write(title + "\n")
        f.write("%5d\n" % num)
        count = 0
        mon_count = 0

        for i in range(len(self.backbone_atoms)):
                resnumber = i + 1
                resname = self.backbone_atoms[i][3]
                tipe = self.backbone_atoms[i][4]
                s = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (
                    resnumber, resname, tipe, resnumber, self.position[self.backbone_atoms[i][0] - 1][0], self.position[self.backbone_atoms[i][0] - 1][1],
                    self.position[self.backbone_atoms[i][0] - 1][2], 0, 0, 0)
                f.write(s)
        if length is None:
            length = 2 * self.max_radius() + 2
        s = "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n" % \
            (length, length, length, 0, 0, 0, 0, 0, 0)
        # (box[X][X], box[Y][Y], box[Z][Z], box[X][Y], box[X][Z], box[Y][X], box[Y][Z], box[Z][X], box[Z][Y])
        f.write(s)
        f.close()

    def go_topology_file(self, title):

        f = open(title + ".top", 'w')
        itp_names = ["go_" + self.itp]
        itp_string = ""
        f2 = open(itp_names[0], 'w')
        f.write('#include "' + "ff.itp" + '"' + "\n")
        for itp_name in itp_names:
            f.write('#include "' + itp_name + '"' + "\n")
        f.write("\n")
        f.write("[ system ]" + "\n")
        f.write(title + "\n")

        #f.write("[ nonbond_params ] \n")
        #f.write("  BB    BB      1   0.86233E-01     0.92953E-03 ; repulsive\n")
        f.write("[ molecules ]\n")
        chain_name = "go_" + self.itp[:-4]
        itp_string += "\n"
        itp_string += self.full_itp_string(chain_name, None, None, None, None, None)
        f.write(chain_name + " " + "1" + "\n")
        f2.write(itp_string)
        f.close()
        f2.close()

    def itp_bonds_string(self, b_obj, num=2):

        string = ""
        for ind, bond in enumerate(self.backbone_bonds + self.backbone_disulfide_bonds):
            string += self.martini_top_string(bond, num=num)
        return string

    def itp_angles_string(self, b_obj, num=3):

        string = ""
        for ind, bond in enumerate(self.backbone_angles):
            string += self.martini_top_string(bond, num=num)
        return string

    def itp_dihedrals_string(self, d_obj, num=4):

        string = ""
        for ind, bond in enumerate(self.backbone_dihedrals):
            string += self.martini_top_string(bond, num=num)
        return string

    def itp_constraints_string(self, c_obj, num=2):
        string = "\n"
        return string

    def itp_restraints_string(self, p_obj, num=1):
        return "\n"

    def itp_pairs_string(self, num=2):

        string = ""
        for ind, bond in enumerate(self.contact_pairs):
            string += self.martini_top_string(bond, num=num)
        return string


    def itp_atoms_string(self, name):

        string = ""
        for ind, atom in enumerate(self.backbone_atoms):
            nr = ind + 1
            tipe = atom[4]
            resnum = ind + 1
            residue = atom[3]
            atomname = tipe
            cnr = ind + 1
            #charge = np.sum(self.charge[self.monomer_indexes[ind]])
            charge = 0
            mass = np.sum([self.mass[x] for x in self.monomer_indexes[ind]])
            string += "%s%5s%5d%8s%5s%5d%8.3f%8.3f\n" % (nr, tipe, resnum, residue, atomname, cnr, charge, mass)
        return string


    def martini_top_string(self, bond, num=2):

        tring = ""
        copied = cp.deepcopy(bond)
        for i in range(num):
            copied[i] = self.backbone_indices.index(bond[i]) + 1
        for i in range(len(bond)):
            tring += str(copied[i]) + " "

        return tring + "\n"

class MartiniProteinFromData(MartiniProtein):

    def __init__(self, positions, types, mass, charge=None, body=None, top=False):

        self.input_position = positions
        self.input_type = types
        self.input_mass = mass
        PDBname = "no_pdb"
        itp = "no_itp"
        super(MartiniProteinFromData, self).__init__(self, PDBname, itp, top=top)

    def build_chain(self, sequence, mma=False, with_ion=False, ps=False):

        self.sequence = []
        self.backbone_atoms = []


        total_count = 0
        for ind in range(len(self.input_position)):
            self.position.append(self.input_position[ind])
            self.type.append(self.input_type[ind])
            self.mass.append(self.input_mass[ind])
            self.charge.append(0)
            self.body.append(-1)
            self.monomer_indexes.append([ind])
            self.sequence.append(self.input_type[ind])
            self.backbone_atoms.append([ind + 1, self.input_type[ind], ind + 1, self.input_type[ind], self.input_type[ind], ind+1, 0])
            total_count += 1
            self.num_beads += 1




