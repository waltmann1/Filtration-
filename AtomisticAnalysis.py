from __future__ import division
import numpy as np
import numpy.linalg as la
import os.path
import networkx as nx
import copy as cp
from matplotlib import pyplot as plt
from matplotlib  import cm
plt.rcParams.update({'font.size': 22})
from mpl_toolkits.mplot3d import Axes3D
import math as m
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib as mpl
import MDAnalysis.transformations as mdat
from Quaternion import QuaternionBetween


class AtomisticAnalysis(object):

    def __init__(self, gro_name, xtc_name, itps = None):

        self.gro_name = gro_name
        self.xtc_name = xtc_name
        self.map = mda.Universe(gro_name)
        self.universe = mda.Universe(xtc_name)
        self.trajectory = self.universe.trajectory

        self.chain_atom_groups = []
        self.trimer_switch = False

        count = 0
        for itp in itps:
            tipes = self.read_itp(itp)
            indices = [int(i + count) for i in range(len(tipes))]
            self.chain_atom_groups.append(self.universe.atoms[indices])
            count += len(tipes)


    def read_itp(self, itp_name):
        f = open(itp_name)
        lines = f.readlines()

        start = False
        tipes = []
        for line in lines:
            s = line.split()
            if len(s) == 0:
                start = False
            elif s[0].isnumeric() and start:
                tipes.append(s[1])
            elif s[:3] == ["[", "atoms", "]"]:
                start = True
        return tipes

    def pentamer_hexamer_angle(self, frame=0):

        time = self.trajectory[frame].time
        #hex_coms = [ self.chain_atom_groups[i].centroid() for i in [5,7,9]]
        #pent_coms = [self.chain_atom_groups[i].centroid() for i in [0, 2, 4]]


        hex_coms = [self.chain_atom_groups[i].centroid() for i in [0, 2, 4]]

        pent_coms = [self.chain_atom_groups[i].centroid() for i in [6, 8, 10]]

        vec1 = np.subtract(hex_coms[1], hex_coms[2])
        vec2 = np.subtract(hex_coms[1], hex_coms[0])

        hex_normal = np.cross(vec1, vec2)



        vec1 = np.subtract(pent_coms[1], pent_coms[2])
        vec2 = np.subtract(pent_coms[1], pent_coms[0])

        pent_normal = np.cross(vec1, vec2)
        #print(pent_coms)

        upn = np.divide(pent_normal, np.linalg.norm(pent_normal))
        uhn = np.divide(hex_normal, np.linalg.norm(hex_normal))
        #print(uhn, upn)

        total_angle = np.arccos(np.dot(upn, uhn))

        upn_bend = cp.deepcopy(upn)
        upn_twist = cp.deepcopy(upn)
        upn_bend[1] = 0
        upn_bend = np.divide(upn_bend, np.linalg.norm(upn_bend))
        upn_twist[0] = 0
        upn_twist = np.divide(upn_twist, np.linalg.norm(upn_twist))


        uhn_bend = cp.deepcopy(uhn)
        uhn_twist = cp.deepcopy(uhn)
        uhn_bend[1] = 0
        uhn_bend = np.divide(uhn_bend, np.linalg.norm(uhn_bend))
        uhn_twist[0] = 0
        uhn_twist = np.divide(uhn_twist, np.linalg.norm(uhn_twist))

        bend_angle = np.arccos(np.dot(upn_bend, uhn_bend))
        twist_angle = np.arccos(np.dot(upn_twist, uhn_twist))

        pointing_vector = np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))

        u_pointing_vector = np.divide(pointing_vector, np.linalg.norm(pointing_vector))

        point_par_vector = np.multiply(np.dot(u_pointing_vector, upn), upn)

        point_perp_vector = np.subtract(u_pointing_vector, point_par_vector)

        u_pointing_perp_vector = np.divide(point_perp_vector, np.linalg.norm(point_perp_vector))

        #print("ahh")
        third_vector = np.cross(upn, u_pointing_perp_vector)

        total_angle = np.arccos(np.dot(upn, uhn))
        bend_angle = np.arccos(np.dot(u_pointing_perp_vector, uhn)) - np.pi / 2
        twist_angle = np.pi / 2 - np.arccos(np.dot(third_vector, uhn))

        #print(uhn_bend, upn_bend)
        #print(uhn_twist, upn_twist)
        return total_angle, bend_angle, twist_angle

    def hexamer_hexamer_angle(self, frame=0):

        time = self.trajectory[frame].time
        #hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in [(0,1), (2,3), (4,5)]]
        hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                    [(6, 7), (8, 9), (10, 11)]]

        vec1 = np.subtract(hex_coms[1], hex_coms[2])
        vec2 = np.subtract(hex_coms[1], hex_coms[0])

        hex_normal = np.cross(vec1, vec2)

        #pent_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in [(6,7), (8,9), (10,11)]]
        pent_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                     [(0, 1), (2, 3), (4, 5)]]

        vec1 = np.subtract(pent_coms[1], pent_coms[2])
        vec2 = np.subtract(pent_coms[1], pent_coms[0])

        pent_normal = np.cross(vec1, vec2)

        upn = np.divide(pent_normal, np.linalg.norm(pent_normal))
        uhn = np.divide(hex_normal, np.linalg.norm(hex_normal))
        #print(uhn, upn)
        #quit()

        total_angle = np.arccos(np.dot(upn, uhn))

        upn_bend = cp.deepcopy(upn)
        upn_twist = cp.deepcopy(upn)
        upn_bend[0] = 0
        upn_bend = np.divide(upn_bend, np.linalg.norm(upn_bend))
        upn_twist[1] = 0
        upn_twist = np.divide(upn_twist, np.linalg.norm(upn_twist))

        uhn_bend = cp.deepcopy(uhn)
        uhn_twist = cp.deepcopy(uhn)
        uhn_bend[0] = 0
        uhn_bend = np.divide(uhn_bend, np.linalg.norm(uhn_bend))
        uhn_twist[1] = 0
        uhn_twist = np.divide(uhn_twist, np.linalg.norm(uhn_twist))

        bend_angle = np.arccos(np.dot(upn_bend, uhn_bend))
        #if upn_bend[0]* uhn_bend[0] < 0:
        #    bend_angle = -bend_angle
        twist_angle = np.arccos(np.dot(upn_twist, uhn_twist))

        # print(uhn_bend, upn_bend)
        # print(uhn_twist, upn_twist)

        pointing_vector = np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))

        u_pointing_vector = np.divide(pointing_vector, np.linalg.norm(pointing_vector))

        point_par_vector = np.multiply(np.dot(u_pointing_vector, upn), upn)

        point_perp_vector = np.subtract(u_pointing_vector, point_par_vector)

        u_pointing_perp_vector = np.divide(point_perp_vector, np.linalg.norm(point_perp_vector))

        print("ahh")
        third_vector = np.cross(upn, u_pointing_perp_vector)

        total_angle = np.arccos(np.dot(upn, uhn))
        bend_angle = np.arccos(np.dot(u_pointing_perp_vector, uhn)) - np.pi / 2
        twist_angle = np.pi / 2 - np.arccos(np.dot(third_vector, uhn))


        return total_angle, bend_angle, twist_angle

    def trimer_hexamer_angle(self, frame=0):

        time = self.trajectory[frame].time
        #hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
        #            [(3, 4), (5, 6), (7, 8)]]

        hex_coms = [self.chain_atom_groups[i].centroid() for i in [3, 6, 7]]

        #hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
         #                      [(0, 1), (2, 3), (4, 5)]]

        vec1 = np.subtract(hex_coms[1], hex_coms[2])
        vec2 = np.subtract(hex_coms[1], hex_coms[0])

        hex_normal = np.cross(vec1, vec2)

        pent_coms = [self.chain_atom_groups[i].centroid() for i in [0,1, 2]]
        #pent_coms = [self.chain_atom_groups[i].centroid() for i in [6, 7, 8]]

        vec1 = np.subtract(pent_coms[1], pent_coms[2])
        vec2 = np.subtract(pent_coms[1], pent_coms[0])

        pent_normal = np.cross(vec1, vec2)
        #print(hex_coms, pent_coms)

        upn = np.divide(pent_normal, np.linalg.norm(pent_normal))
        uhn = np.divide(hex_normal, np.linalg.norm(hex_normal))
        #print(uhn, upn)

        #print(pent_coms, hex_coms)
        pointing_vector = np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))

        u_pointing_vector = np.divide(pointing_vector, np.linalg.norm(pointing_vector))

        point_par_vector = np.multiply(np.dot(u_pointing_vector, upn), upn)

        point_perp_vector = np.subtract(u_pointing_vector, point_par_vector)

        u_pointing_perp_vector = np.divide(point_perp_vector, np.linalg.norm(point_perp_vector))

        third_vector = np.cross(upn, u_pointing_perp_vector)

        total_angle = np.arccos(np.dot(upn, uhn))

        """
        upn_bend = cp.deepcopy(upn)
        upn_twist = cp.deepcopy(upn)
        upn_bend[1] = 0
        upn_bend = np.divide(upn_bend, np.linalg.norm(upn_bend))
        upn_twist[0] = 0
        upn_twist = np.divide(upn_twist, np.linalg.norm(upn_twist))


        uhn_bend = cp.deepcopy(uhn)
        uhn_twist = cp.deepcopy(uhn)
        uhn_bend[1] = 0
        uhn_bend = np.divide(uhn_bend, np.linalg.norm(uhn_bend))
        uhn_twist[0] = 0
        uhn_twist = np.divide(uhn_twist, np.linalg.norm(uhn_twist))
        """
        #bend_angle = np.arccos(np.dot(upn_bend, uhn_bend))
        #twist_angle = np.arccos(np.dot(upn_twist, uhn_twist))

        bend_angle = np.arccos(np.dot(u_pointing_perp_vector, uhn)) - np.pi/2
        twist_angle = np.pi/2 - np.arccos(np.dot(third_vector, uhn))

        #print(uhn_bend, upn_bend)
        #print(uhn_twist, upn_twist)
        return total_angle, bend_angle, twist_angle


    def t_trimer_hexamer_angle(self, frame=0):

        if not self.trimer_switch:
            temp = cp.deepcopy(self.chain_atom_groups)

            new = temp[3:]
            new.extend(temp[:3])
            self.chain_atom_groups = new
            self.trimer_switch = True
        return self.trimer_hexamer_angle(frame=frame)


    def z_dist(self, frame=0):

        time = self.trajectory[frame].time
        hex_coms = [self.chain_atom_groups[i].centroid() for i in [1, 3, 5]]

        pent_coms = [self.chain_atom_groups[i].centroid() for i in [6, 8, 10]]

        return np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))[2]

    def z_dist_pent(self, frame=0):

        time = self.trajectory[frame].time
        #hex_coms = [self.chain_atom_groups[i].centroid() for i in [9, 7, 5]]

        hex_coms = [self.chain_atom_groups[i].centroid() for i in [0, 2, 4]]


        pent_coms = [self.chain_atom_groups[i].centroid() for i in [10, 8, 6]]
        #pent_coms = [self.chain_atom_groups[i].centroid() for i in [0, 2, 4]]


        return np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))[2]


    def y_dist(self, frame=0):

        time = self.trajectory[frame].time
        hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                    [(0, 1), (2, 3), (4, 5)]]
        pent_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                     [(6, 7), (8, 9), (10, 11)]]

        return np.subtract(np.average(hex_coms, axis=0), np.average(pent_coms, axis=0))[1]


    def com_dist(self, frame=0):

        time = self.trajectory[frame].time
        hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                    [(0, 1), (2, 3), (4, 5)]]
        pent_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                     [(6, 7), (8, 9), (10, 11)]]

        return np.linalg.norm(np.subtract(np.average(hex_coms, axis=0), np.average(pent_coms, axis=0)))



    def z_dist(self, frame=0):

        time = self.trajectory[frame].time
        hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                    [(0, 1), (2, 3), (4, 5)]]
        pent_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
                     [(6, 7), (8, 9), (10, 11)]]

        return np.subtract(np.average(hex_coms, axis=0), np.average(pent_coms, axis=0))[2]


    def z_dist_trimer(self, frame=0):

        time = self.trajectory[frame].time
        #hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
         #           [(0, 1), (2, 3), (4, 5)]]


        #pent_coms = [self.chain_atom_groups[i].centroid() for i in [6, 7, 8]]

        #hex_coms = [(self.chain_atom_groups[i[0]] + self.chain_atom_groups[i[1]]).centroid() for i in
        #          [(3, 4), (5, 6), (7, 8)]]

        hex_coms = [self.chain_atom_groups[i].centroid() for i in [3, 4, 5, 6, 7, 8]]

        pent_coms = [self.chain_atom_groups[i].centroid() for i in [0, 1, 2]]

        return np.subtract(np.average(pent_coms, axis=0), np.average(hex_coms, axis=0))[2]

    def vector_from_indices(self, index1, index2, frame=0):


        time = self.trajectory[frame].time
        atom1 = self.universe.select_atoms("index " + str(index1))[0]
        atom2 = self.universe.select_atoms("index " + str(index2))[0]
        pos1 = atom1.position
        pos2 = atom2.position
        vec = np.subtract(pos2, pos1)
        #print(atom1.position, atom2.position, vec)
        return vec

    def lys_angle(self, c1index, n1index, c2index, n2index, frame=0):

        vec1 = self.vector_from_indices(c1index, n1index, frame=frame)
        vec2 = self.vector_from_indices(c2index, n2index, frame=frame)
        u1 = np.divide(vec1, np.linalg.norm(vec1))
        u2 = np.divide(vec2, np.linalg.norm(vec2))
        angle = np.arccos(np.dot(u1,u2))
        return np.rad2deg(angle)

