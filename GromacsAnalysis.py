from __future__ import division
import numpy as np
import numpy.linalg as la
import os.path
#import networkx as nx
import copy as cp
from matplotlib import pyplot as plt
from matplotlib  import cm
import matplotlib
plt.rcParams.update({'font.size': 22})
from mpl_toolkits.mplot3d import Axes3D
import math as m
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib as mpl
import MDAnalysis.transformations as mdat
from Quaternion import QuaternionBetween
from scipy.optimize import curve_fit
import scipy.special as sc
import math

class GromacsAnalysis(object):

    def __init__(self, gro_name, xtc_name, unwrap=False, count=1, itps = None):

        self.gro_name = gro_name
        self.xtc_name = xtc_name
        self.map = mda.Universe(gro_name)
        self.universe = mda.Universe(xtc_name)
        self.trajectory = self.universe.trajectory
        if unwrap:
            ag = self.universe.atoms
            transform = mdat.unwrap(ag)
            self.trajectory.add_transformations(transform)
        self.chains = []
        self.map_chains = []
        #count = 1
        map_atoms = self.map.select_atoms("resname Po" + str(count))
        subcount = 0
        self.num_protein_chains = 0
        if len(map_atoms) == 0:
            self.num_protein_chains = 1
            try:
                while self.map.atoms[subcount].resname[:2] != "Po":
                    subcount += 1
                map_atoms = self.map.atoms[:subcount]
            except IndexError:
                print("No polymers")


        indices = map_atoms.indices
        atoms = self.universe.atoms[indices]
        while len(atoms) > 0:
            self.chains.append(atoms)
            self.map_chains.append(map_atoms)
            count += 1
            map_atoms = self.map.select_atoms("resname Po" + str(count))
            indices = map_atoms.indices
            atoms = self.universe.atoms[indices]
        #print(self.chains, len(self.chains))
        if count > 1 and itps is not None:
            new = []
            map_new = []
            lengths = [len(self.read_itp(itp)) for itp in itps]
            start = 0
            for length in lengths:
                map_atoms = self.map.atoms[start:start + length]
                indices = map_atoms.indices
                atoms = self.universe.atoms[indices]
                new.append(atoms)
                map_new.append(map_atoms)
                start += length
            new.extend(self.chains[1:])
            map_new.extend(self.map_chains[1:])
            self.chains = new
            self.map_chains = map_new
            self.num_protein_chains = len(itps)

        self.chains_monomer_list = None
        self.chains_beads_by_monomer = None



    def monomerize_mma_chain(self, chain_atom_group):

        chain_by_monomers = []
        atoms_with_monomer = []
        index = 0
        while chain_atom_group[index].name == "SC1":
            chain_by_monomers.append([index])
            atoms_with_monomer.append([index, 0, -1])
            index += 1
        #print("a", chain_by_monomers)
        num_monomers = len(chain_by_monomers)
        mon_index = -1
        index_in_monomer = 1
        while index < len(chain_atom_group.indices):
            if chain_atom_group[index].name == "Na":
                mon_index += 1
                index_in_monomer = 1
            #print(mon_index, index)
            chain_by_monomers[mon_index].append(index)
            atoms_with_monomer.append([mon_index, index_in_monomer, -1])
            index += 1
            index_in_monomer += 1

        names = []
        for monomer in range(num_monomers):
            if len(chain_by_monomers[monomer]) == 11:
                names.append("PEOPEMA")
            if len(chain_by_monomers[monomer]) == 4:
                atom = chain_atom_group[chain_by_monomers[monomer][3]]
                if atom.name == "Qa":
                    names.append("SPMA")
                else:
                    names.append("EHMA")
        for bead in range(len(atoms_with_monomer)):
            atoms_with_monomer[bead][2] = names[atoms_with_monomer[bead][0]]

        return names, atoms_with_monomer


    def radius_of_gyration(self, index, frame=0):
        group = self.chains[index]
        return self.group_radius_of_gyration(group, frame=frame)

    def chain_rdf(self, index, frames):

        names = [self.map_chains[index][i].name for i in range(len(self.map_chains[index]))]

        unique_names = list(np.unique(names))

        rdf_distances = [[] for _ in range(len(unique_names))]

        for frame in frames:
            group = self.chains[index]
            #positions = np.array(self.get_positions_chain(index, frame))
            positions = group.positions
            #com = np.average(positions)
            com = self.group_center_of_mass(group, frame)
            #dist = la.norm(np.subtract(com, positions), axis=1)
            dist = distances.distance_array(positions, com, box=self.universe.dimensions)
            for ind, d in enumerate(dist):
                the_index = unique_names.index(names[ind])
                rdf_distances[the_index].append(d)

        return unique_names, rdf_distances

    def rdf_polymers_petase_center(self, indices, frames):

        names = [self.map_chains[index][i].name for index in indices for i in range(len(self.map_chains[index]))]

        unique_names = list(np.unique(names))

        rdf_distances = [[] for _ in range(len(unique_names))]

        group = self.chains[indices[0]]
        for index in indices[1:]:
            group = group.concatenate(self.chains[index])

        for frame in frames:
            # positions = np.array(self.get_positions_chain(index, frame))
            positions = group.positions
            # com = np.average(positions)
            com = self.chains[0].positions[264]
            # dist = la.norm(np.subtract(com, positions), axis=1)
            dist = distances.distance_array(positions, com, box=self.universe.dimensions)
            for ind, d in enumerate(dist):
                the_index = unique_names.index(names[ind])
                rdf_distances[the_index].append(d)

        return unique_names, rdf_distances

    def chains_rdf(self, indices, frames):

        names = [self.map_chains[index][i].name for index in indices for i in range(len(self.map_chains[index]))]

        unique_names = list(np.unique(names))

        rdf_distances = [[] for _ in range(len(unique_names))]

        group = self.chains[indices[0]]

        for chain in self.chains[1:]:
            group = group + chain

        for frame in frames:
            #positions = np.array(self.get_positions_chain(index, frame))
            positions = group.positions
            #com = np.average(positions)
            com = self.group_center_of_mass(group, frame)
            #dist = la.norm(np.subtract(com, positions), axis=1)
            dist = distances.distance_array(positions, com, box=self.universe.dimensions)
            for ind, d in enumerate(dist):
                the_index = unique_names.index(names[ind])
                rdf_distances[the_index].append(d)

        return unique_names, rdf_distances

    def z_dist_chain(self, index, frame):


        group = self.chains[index]
        time = self.trajectory[frame].time
        z_positions = group.positions[:,2]

        return z_positions

    def z_dist_chain_types(self, index, frame):

        names = [self.map_chains[index][i].name for i in range(len(self.map_chains[index]))]
        unique_names = sorted(list(np.unique(names)))

        group = self.chains[index]
        time = self.trajectory[frame].time
        z_positions = group.positions[:, 2]
        rdf_distances = [[] for _ in range(len(unique_names))]

        for ind, d in enumerate(z_positions):
            the_index = unique_names.index(names[ind])
            rdf_distances[the_index].append(d)

        return unique_names, rdf_distances


    def graph_hist_z_chaintypes(self, index, frames, norm=True, save_name="z_hist_types.png"):


        u_names, rdf_distances = self.z_dist_chain_types(index, frames[0])

        print(u_names)
        print(rdf_distances[0])

        distances = list(self.z_dist_chain(index, frames[0]))

        for frame in frames[1:]:
            new_u_names, new_rdf_distances = self.z_dist_chain_types(index, frame)
            for ind in range(len(rdf_distances)):
                rdf_distances[ind].extend(new_rdf_distances[ind])


        bins = [1 * i for i in range(400)]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('z (nm)')
        ax1.set_ylabel('P(z)')



        for ind in range(len(rdf_distances)):
            hist, bin_edges = np.histogram(rdf_distances[ind], bins=bins)
            if norm:
                hist = np.divide(hist, np.sum(hist))
            bin_edges = np.divide(bin_edges, 10)
            ax1.plot(bin_edges[:-1], hist, label= u_names[ind])

        plt.legend(fontsize=14, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def name_z_pos(self, name, frame):

        map_atoms = self.map.select_atoms("name " + name)
        indices = map_atoms.indices
        atoms = self.universe.atoms[indices]

        time = self.trajectory[frame].time
        return atoms.positions[:, 2]


    def graph_interface(self, frame, water_name="WP", oil_name="CX", save_name="interface_hist.png"):

        oil_pos = self.name_z_pos(oil_name, frame)
        water_pos = self.name_z_pos(water_name, frame)

        bins = [1 * i for i in range(400)]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_title("Frame " + str(frame))
        ax1.set_xlabel('z (nm)')
        ax1.set_ylabel('P(z)')

        oil_hist, bin_edges = np.histogram(oil_pos, bins=bins)
        oil_hist = np.divide(oil_hist, np.sum(oil_hist))
        bin_edges = np.divide(bin_edges, 10)
        ax1.plot(bin_edges[:-1], oil_hist, label=oil_name)

        water_hist, bin_edges = np.histogram(water_pos, bins=bins)
        water_hist = np.divide(water_hist, np.sum(water_hist))
        bin_edges = np.divide(bin_edges, 10)
        ax1.plot(bin_edges[:-1], water_hist, label=water_name)



        plt.legend(fontsize=14, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()



    def graph_hist_z(self, indices, frames, save_name="z_hist.png"):

        distances = [[] for i in range(len(indices))]

        for ind, chain_index in enumerate(indices):
            distances[ind] = list(self.z_dist_chain(chain_index, frames[0]))

            for frame in frames[1:]:
                distances[ind].extend(list(self.z_dist_chain(chain_index, frame)))


        bins = [1 * i for i in range(400)]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('z (nm)')
        ax1.set_ylabel('P(z)')



        for ind in range(len(indices)):
            hist, bin_edges = np.histogram(distances[ind], bins=bins)
            box_area = 400
            z_space = bin_edges[1] - bin_edges[0]
            num_frames = len(frames)
            total_factor = box_area * z_space * num_frames
            hist = np.divide(hist, total_factor)
            #hist = np.divide(hist, np.sum(hist))
            bin_edges = np.divide(bin_edges, 10)
            ax1.plot(bin_edges[:-1], hist, label="Chain Index: " + str(indices[ind]))

        plt.legend(fontsize=14, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_petase_hist(self, index, frames):


        unique_names, rdf_distances = self.rdf_polymers_petase_center(index, frames)

        bins = [1 * i for i in range(0, 200, 2)]
        fig = plt.figure()

        #print(unique_names)
        #print(distances)

        #quit()
        #combined = cp.deepcopy(rdf_distances[0])
        #combined.extend(rdf_distances[2])
        #combined.extend(rdf_distances[4])

        #unique_names.append("combined")
        #rdf_distances.append(combined)

        #combined2 = cp.deepcopy(rdf_distances[1])
        #combined2.extend(rdf_distances[5])

        #unique_names.append("combined2")
        #rdf_distances.append(combined2)

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('R (nm)')
        ax1.set_ylabel('P(R)')

        myorder = [6, 7, 3, 5, 1, 0, 4, 2]

        np_list = [0, 2, 4, 6]

        names = ["NonPolar (C1)", "Polar (EO)", "NonPolar (Na)", "Negatively Charged (Qa)", "NonPolar (SC1)", "Polar (SP2)", "NonPolar (All)", "Polar (All)"]
        colors = ["#aaaa00",  "#0000ff", "#aaaa00", "#00ff7f", "#aaaa00",  "#0000ff", "#aaaa00", "#0000ff"]
        shapes = ["dashed", "dashed", "dotted", "solid", "dashdot", "dotted", "solid", "solid"]

        #for ind in myorder:
        for ind in range(len(unique_names)):
            hist, bin_edges = np.histogram(rdf_distances[ind], bins=bins)
            #if ind in np_list:
                #hist = np.divide(hist, 2)
            hist = np.divide(hist, np.sum(hist))
            bin_edges = np.divide(bin_edges, 10)
            #if shapes[ind] == "solid":
             #   mpl.rcParams['lines.linewidth'] = 2
            #else:
            #    mpl.rcParams['lines.linewidth'] = 1
            #ax1.plot(bin_edges[:-1], hist, label=names[ind], color=colors[ind], linestyle=shapes[ind])
            ax1.plot(bin_edges[:-1], hist, label=unique_names[ind], color=colors[ind], linestyle=shapes[ind])

        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig("petase_hist.png", bbox_inches='tight', pad_inches=.2)
        plt.show()


    def graph_hist(self, index, frames):


        unique_names, rdf_distances = self.chain_rdf(index, frames)
        bins = [1 * i for i in range(0,52, 2)]
        fig = plt.figure()

        #print(unique_names)
        #print(distances)

        #quit()
        combined = cp.deepcopy(rdf_distances[0])
        combined.extend(rdf_distances[2])
        combined.extend(rdf_distances[4])

        unique_names.append("combined")
        rdf_distances.append(combined)

        combined2 = cp.deepcopy(rdf_distances[1])
        combined2.extend(rdf_distances[5])

        unique_names.append("combined2")
        rdf_distances.append(combined2)

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('R (nm)')
        ax1.set_ylabel('P(R)')

        myorder = [6, 7, 3, 5, 1, 0, 4, 2]

        np_list = [0, 2, 4, 6]

        names = ["NonPolar (C1)", "Polar (EO)", "NonPolar (Na)", "Negatively Charged (Qa)", "NonPolar (SC1)", "Polar (SP2)", "NonPolar (All)", "Polar (All)"]
        colors = ["#aaaa00",  "#0000ff", "#aaaa00", "#00ff7f", "#aaaa00",  "#0000ff", "#aaaa00", "#0000ff"]
        shapes = ["dashed", "dashed", "dotted", "solid", "dashdot", "dotted", "solid", "solid"]

        for ind in myorder:
            hist, bin_edges = np.histogram(rdf_distances[ind], bins=bins)
            #if ind in np_list:
                #hist = np.divide(hist, 2)
            hist = np.divide(hist, np.sum(hist))
            bin_edges = np.divide(bin_edges, 10)
            if shapes[ind] == "solid":
                mpl.rcParams['lines.linewidth'] = 2
            else:
                mpl.rcParams['lines.linewidth'] = 1
            ax1.plot(bin_edges[:-1], hist, label=names[ind], color=colors[ind], linestyle=shapes[ind])

        #plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig("micelle_hist.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def group_diameter(self, index, frame=0):

        time = self.trajectory[frame].time
        group = self.chains[index]
        #positions = np.array(self.get_positions_chain(index, frame))
        positions = group.positions
        self_distances = distances.self_distance_array(positions, box=self.universe.dimensions)

        return np.max(self_distances)

    def group_radius_of_gyration(self, group, frame=0):

        time = self.trajectory[frame].time
        self_distances = distances.self_distance_array(group.positions, box=self.universe.dimensions)
        sq = np.multiply(self_distances, self_distances)
        ave = np.average(sq, axis=0)
        msrg = ave / 6
        rg = np.sqrt(msrg)
        return rg

    def group_center_of_mass(self, group, frame=0):

        time = self.trajectory[frame].time

        return  group.center_of_geometry(pbc=True)
        #mass_array = self.map.atoms[group.indices].masses
        #position_array = group.positions
        #array = [0, 0, 0]

        #for i in range(len(position_array)):
        #    for x in range(3):
        #        array[x] += mass_array[i] * position_array[i][x]

        #return array / np.sum(mass_array)

    def positions_gyration_tensor(self, positions):

        gyr_ten = np.zeros((3, 3))
        # print("positions", positions)
        com = np.average(positions, axis=1)
        # print("com", com)

        for pos in positions:
            gyr_ten[0][0] += (pos[0] - com[0]) ** 2
            gyr_ten[1][0] += (pos[0] - com[0]) * (pos[1] - com[1])
            gyr_ten[2][0] += (pos[0] - com[0]) * (pos[2] - com[2])
            gyr_ten[0][1] += (pos[1] - com[1]) * (pos[0] - com[0])
            gyr_ten[1][1] += (pos[1] - com[1]) ** 2
            gyr_ten[2][1] += (pos[1] - com[1]) * (pos[2] - com[2])
            gyr_ten[0][2] += (pos[2] - com[2]) * (pos[0] - com[0])
            gyr_ten[1][2] += (pos[2] - com[2]) * (pos[1] - com[1])
            gyr_ten[2][2] += (pos[2] - com[2]) ** 2
        return gyr_ten


    def group_gyration_tensor(self, group, frame=0):

        positions = group.positions
        gyr_ten = np.zeros((3,3))
        #print("positions", positions)
        com = self.group_center_of_mass(group, frame=frame)
        #print("com", com, len(positions))

        for pos in positions:
            gyr_ten[0][0] += (pos[0] - com[0]) ** 2
            gyr_ten[1][0] += (pos[0] - com[0]) * (pos[1] - com[1])
            gyr_ten[2][0] += (pos[0] - com[0]) * (pos[2] - com[2])
            gyr_ten[0][1] += (pos[1] - com[1]) * (pos[0] - com[0])
            gyr_ten[1][1] += (pos[1] - com[1]) ** 2
            gyr_ten[2][1] += (pos[1] - com[1]) * (pos[2] - com[2])
            gyr_ten[0][2] += (pos[2] - com[2]) * (pos[0] - com[0])
            gyr_ten[1][2] += (pos[2] - com[2]) * (pos[1] - com[1])
            gyr_ten[2][2] += (pos[2] - com[2]) ** 2
        return gyr_ten

    def group_asphericity(self, group, frame=0):

        gyr_ten = self.group_gyration_tensor(group, frame=frame)
        #print("gyr", gyr_ten)
        eigs = np.linalg.eigvals(gyr_ten)
        eigs.sort()
        print("eigs", eigs)
        print(eigs[2]- eigs[1] , (la.norm(eigs) **2))
        return (eigs[2] ** 2 - .5 * (eigs[0] ** 2 + eigs[1] ** 2)) / (la.norm(eigs) **2)

    def group_acylindricity(self, group, frame=0):

        gyr_ten = self.group_gyration_tensor(group, frame=frame)
        eigs = np.linalg.eigvals(gyr_ten)
        eigs.sort()
        print("acyl eigs", eigs)
        print((eigs[2] ** 2 - .5 * (eigs[0] ** 2 + eigs[1] ** 2)) / (la.norm(eigs) **2) )
        return (eigs[1] ** 2 - eigs[0]**2) / (eigs[1]**2 + eigs[0]**2)

    def check_contact(self, index_1, index_2, cut=5.3, frame=0):

        time = self.trajectory[frame].time
        dist_arr = distances.distance_array(self.chains[index_1].positions, self.chains[index_2].positions, box=self.universe.dimensions)
        if np.amin(np.array(dist_arr)) < cut:
            return True
        return False

    def get_all_contacts_multiples(self, indices_1, indices_2, cut=5.3, frame=0):

        for index_1 in indices_1:
            for index_2 in indices_2:
                self.get_all_contacts(index_1, index_2, cut=cut, frame=frame)

    def average_pet_contact_groups(self, frames, pet_indices):

        types = []
        polar_count = 0
        charged_count = 0
        non_polar_count = 0
        for frame in frames:
            print(frame)
            for ind, chain in enumerate(self.chains):
                if ind not in pet_indices:
                    for pet_chain_index in pet_indices:
                        con = self.get_all_contacts(pet_chain_index, ind, frame=frame)
                        types.extend(self.contacts_types(con, pet_chain_index, ind))
        typed, counted = self.count_contacts(types)
        #for i in range(len(typed)):
            #print(typed[i], counted[i])
        for ind in range(len(typed)):
            print(typed[ind], counted[ind])
            if "SC1" in typed[ind][1]:
                polar_count += counted[ind]
            elif "Na" in typed[ind][1]:
                charged_count += counted[ind]
            else:
                non_polar_count += counted[ind]

        return non_polar_count/len(frames), polar_count/len(frames), charged_count/len(frames)



    def get_all_contacts(self, index_1, index_2, cut=5.3, frame=0):

        #print("lets get the contacts")
        if os.path.exists(str(self.xtc_name[:-4]) + "_frame" + str(frame) + "_contacts_" + str(index_1) + "_" + str(index_2) + ".txt"):
            return self.read_data(str(self.xtc_name[:-4]) + "_frame" + str(frame) + "_contacts_"  + str(index_1) + "_" + str(index_2) + ".txt")


        time = self.trajectory[frame].time
        dist_arr = distances.distance_array(self.chains[index_1].positions, self.chains[index_2].positions,
                                            box=self.universe.dimensions)
        pairs = []
        for i in range(dist_arr.shape[0]):
            for j in range(dist_arr.shape[1]):
                if dist_arr[i][j] < cut:
                    pairs.append([i, j])
        self.write_pairs(pairs, index_1, index_2, frame)

        return pairs



    def write_pairs(self, pairs, index_1, index_2, frame):

        write_name = str(self.xtc_name[:-4]) + "_frame" + str(frame) + "_contacts_" + str(index_1) + "_" + str(index_2) + ".txt"

        f = open(write_name, 'w')

        for pair in pairs:
            f.write(str(pair[0]) + " " + str(pair[1]) + "\n")
        f.close()


    def read_data(self, path):

        f = open(path)
        data = f.readlines()
        out = []
        for line in data:
            s = line.split()
            part = [int(thing) for thing in s]
            out.append(part)
        return out

    def protein_times_contacted(self, protein_index, frame=0, cut=5.3, chain_indices=None):

        count = np.zeros(len(self.chains[protein_index].positions))

        if chain_indices is not None:
            for index in chain_indices:
                pairs = self.get_all_contacts(protein_index, index, cut=cut, frame=frame)
                for pair in pairs:
                    count[pair[0]] += 1
            return count


        for index in list(range(self.num_protein_chains, len(self.chains))):
            pairs = self.get_all_contacts(protein_index, index, cut=cut, frame=frame)
            for pair in pairs:
                count[pair[0]] += 1

        return count


    def average_protein_times_contacted(self, frames, cut=5.3, protein_index=0, chain_indices=None):

        count = self.protein_times_contacted(protein_index, frames[0], cut=cut, chain_indices=chain_indices)

        for frame in frames[1:]:
            count = np.add(count, self.protein_times_contacted(protein_index, frame=frame, cut=cut, chain_indices=chain_indices))

        return np.divide(count, len(frames))


    def graph_contacts(self, frames, protein_index=0, cut=5.3, save_name="contact_numbers.png", petase=True, p450=False, chain_indices=None):

            cvs = self.average_protein_times_contacted(frames, cut=cut, protein_index=protein_index, chain_indices=chain_indices)

            fig = plt.figure()

            ax1 = fig.add_subplot(111)

            ax1.set_title("Average Number of Contacts")

            active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                           423, 424, 425, 426, 427, 428]

            f_g_loop = list(range(399,456, 1))


            ax1.set_xlabel('Bead index')
            ax1.set_ylabel('Average number of contacts')
            ax1.plot(list(range(len(cvs))), cvs)
            if petase:
                for act in active_site:
                    ax1.plot([act - 1 for _ in range(10)], [i * .5 for i in range(10)],color="red", label="Active Site")
            if p450:
                for act in f_g_loop:
                    ax1.plot([act - 1 for _ in range(10)], [i * .5 for i in range(10)], color=(1, 0, 0, 0.2),
                             label="F-G Channel")
            #plt.legend()
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
            plt.show()

    def graph_contacts_from_input(self, input, petase=True, p450=False, save_name="contact_numbers_from_input.png", write=False):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")

        active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                       423, 424, 425, 426, 427, 428]

        f_g_loop = list(range(399, 456, 1))

        ax1.set_xlabel('Bead Index')
        ax1.set_ylabel('Average Number of Contacts')
        ax1.plot(list(range(len(input))), input)

        threshold = .95
        #ax1.plot(list(range(len(input))), [threshold for _ in range(len(input))], color="blue", label="Threshold", linestyle="dashed")


        if petase:
            for act in active_site:
                ax1.plot([act for _ in range(10)], [i * .3 for i in range(10)], color="gray",linestyle="dotted", label="Active Site")
        if p450:
            for act in f_g_loop:
                ax1.plot([act - 1 for _ in range(10)], [i * .5 for i in range(10)], color="red", label="F-G Channel")

        if write:
            f = open("protein_index_contact_numbers.txt", 'w')
            for i in range(len(input)):
                f.write(str(i) + " " + str(input[i]) + "\n")
        #plt.legend()
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_spatial_contacts_from_input(self, input, protein_index=0, cut=5.3, save_name="contact_spatial_from_input.png", petase=True, p450=False):

        positions = self.map_chains[protein_index].positions
        com = np.average(positions, axis=0)
        positions = np.subtract(positions, com)
        quat = QuaternionBetween(positions[118], [1, 0, 0])
        for ind, position in enumerate(positions):
            positions[ind] = quat.orient(position)

        w_spherical = self.appendSpherical_np(positions)
        phi_o = w_spherical[:, 4]
        theta_o = w_spherical[:, 5]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_title("Spatial Contact Distribution")

        ax1.set_xlabel('$\u03B8$', fontsize=30)
        ax1.set_ylabel('$\phi$', fontsize=30)
        to_remove = []
        greater = []
        less = []
        for ind in range(len(input)):
            if np.abs(input[ind]) < .95:
                to_remove.append(ind)
            else:
                if input[ind] > .95:
                    greater.append(ind)
                elif input[ind] < -.95:
                    less.append(ind)
        phi = np.delete(phi_o, to_remove)
        theta = np.delete(theta_o, to_remove)
        cvs = np.delete(input, to_remove)
        print(len(phi))

        active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                        423, 424, 425, 426, 427, 428]

        f_g_loop = list(range(399, 456,1))

        ax1.set_xlim(-3.14, 3.14)
        ax1.set_ylim(-np.pi/2, np.pi/2)

        colormap = matplotlib.colors.ListedColormap(["purple", "yellow"])
        ax1.scatter(theta, phi - np.pi/2, c=cvs, cmap=colormap, vmin=.5, vmax=.6, s=100)
        #ax1.scatter(theta, phi - np.pi / 2, c=cvs, cmap=cm.viridis, vmin=.5, vmax=.6, s=100)
        if petase:
            for act in active_site:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act] - np.pi/2, c="gray", s=100)
                else:
                    ax1.scatter(theta_o[act], phi_o[act] - np.pi/2, c="cyan", s=100)

        if p450:
            for act in f_g_loop:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act] , c="red")
                else:
                    ax1.scatter(theta_o[act], phi_o[act] , c="red", marker="x")

        g_string = ""
        l_string = ""
        a_string = ""
        for thing in greater:
            g_string += "index == " + str(thing) + " || "
        for thing in less:
            l_string += "index == " + str(thing) + " || "
        for thing in active_site:
            a_string += "index == " + str(thing) + " || "
        print(g_string)
        print(l_string)
        print(a_string)
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()


    def graph_spatial_contacts(self, frames, protein_index=0, cut=5.3, save_name="contact_spatial.png", petase=True, p450=False, chain_indices=None):

        cvs_o = self.average_protein_times_contacted(frames, cut=cut, protein_index=protein_index, chain_indices=chain_indices)
        positions = self.map_chains[protein_index].positions
        com = np.average(positions, axis=0)
        positions = np.subtract(positions, com)
        quat = QuaternionBetween(positions[118], [1, 0, 0])
        for ind, position in enumerate(positions):
            positions[ind] = quat.orient(position)

        w_spherical = self.appendSpherical_np(positions)
        phi_o = w_spherical[:, 4]
        theta_o = w_spherical[:, 5]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_title("Spatial Contact Distribution")

        ax1.set_xlabel('Theta')
        ax1.set_ylabel('Phi')
        to_remove = []
        for ind in range(len(cvs_o)):
            if cvs_o[ind] < 1:
                to_remove.append(ind)
        phi = np.delete(phi_o, to_remove)
        theta = np.delete(theta_o, to_remove)
        cvs = np.delete(cvs_o, to_remove)
        print(len(phi))

        active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                        423, 424, 425, 426, 427, 428]

        f_g_loop = list(range(399, 456,1))

        ax1.set_xlim(-3.14, 3.14)
        ax1.set_ylim(0, 3.14)
        ax1.scatter(theta, phi, c=cvs, cmap=cm.viridis, vmin=1, vmax=5)
        if petase:
            for act in active_site:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act], c="red")
                else:
                    ax1.scatter(theta_o[act], phi_o[act], c="red", marker="x")

        if p450:
            for act in f_g_loop:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act], c="red")
                else:
                    ax1.scatter(theta_o[act], phi_o[act], c="red", marker="x")
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()


    def graph_spatial_types(self, itp, protein_index=0, save_name="types_spatial.png", petase=True, p450=False):


        positions = self.map_chains[protein_index].positions
        com = np.average(positions, axis=0)
        active_site = [117, 118, 267, 268, 269,270,271,272, 312,313,314,315,316,360,361,
                       423,424,425,426,427,428]

        f_g_loop = list(range(399, 456, 1))

        positions = np.subtract(positions, com)
        if True:
            quat = QuaternionBetween(positions[118], [1, 0, 0])
        else:
            quat = QuaternionBetween(positions[0], [1, 0, 0])
        for ind, position in enumerate(positions):
            positions[ind] = quat.orient(position)
        w_spherical = self.appendSpherical_np(positions)
        phi = w_spherical[:, 4]
        theta = w_spherical[:, 5]
        itp = self.read_itp(itp)
        print(len(itp), "len(itp)" )
        colors = []
        #quit()
        for ind, type in enumerate(itp):
            if ind in active_site and petase:
                colors.append("purple")
            elif ind in f_g_loop and p450:
                colors.append("purple")
            elif "P" in type:
                colors.append("blue")
            elif "Qa" in type:
                colors.append("green")
            elif "Qd" in type:
                colors.append("red")
            elif "C" in type:
                colors.append("k")
            else:
                colors.append("yellow")

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_title("Spatial Type Distribution")

        ax1.set_xlabel('Theta')
        ax1.set_ylabel('Phi')

        print(len(phi))
        ax1.set_xlim(-3.14, 3.14)
        ax1.set_ylim(0, 3.14)
        print(len(theta))
        ax1.scatter(theta, phi, c=colors, label=itp)
        #plt.legend()
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()



    def appendSpherical_np(self, xyz):

        ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
        xy = xyz[:, 0] ** 2 + xyz[:, 1] ** 2
        ptsnew[:, 3] = np.sqrt(xy + xyz[:, 2] ** 2)
        ptsnew[:, 4] = np.arctan2(np.sqrt(xy), xyz[:, 2])  # for elevation angle defined from Z-axis down
        # ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
        ptsnew[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])
        return ptsnew



    def contacts_types(self, pairs, index1, index2, itp1=None, itp2=None):

        type_pairs = []
        if itp1 is not None:
            itp1 = self.read_itp(itp1)
        if itp2 is not None:
            itp2 = self.read_itp(itp2)
        for pair in pairs:
            temp = [self.map_chains[index1][pair[0]].name, self.map_chains[index2][pair[1]].name]
            if itp1 is not None:
                temp[0] = itp1[pair[0]]
            if itp2 is not None:
                temp[1] = itp2[pair[1]]
            type_pairs.append(temp)
        return type_pairs


    def init_mma_monomer_lists(self):

        if self.chains_monomer_list == None:
            self.chains_monomer_list = [[] for _ in range(self.num_protein_chains)]
            self.chains_beads_by_monomer = [[] for _ in range(self.num_protein_chains)]
            for i in range(self.num_protein_chains, len(self.map_chains)):
                names, second = self.monomerize_mma_chain(self.map_chains[i])
                self.chains_monomer_list.append(names)
                self.chains_beads_by_monomer.append(second)


    def contacts_monomer(self, pairs, index1, index2, itp1=None, itp2=None):

        #assuming the second part of the pair is an mma-based monomer

        self.init_mma_monomer_lists()
        type_pairs = []
        if itp1 is not None:
            itp1 = self.read_itp(itp1)
        if itp2 is not None:
            itp2 = self.read_itp(itp2)
        for pair in pairs:
            print(self.chains_beads_by_monomer[index2][pair[1]])
            temp = [self.map_chains[index1][pair[0]].name, self.chains_beads_by_monomer[index2][pair[1]]]
            if itp1 is not None:
                temp[0] = itp1[pair[0]]
            if itp2 is not None:
                temp[1] = itp2[pair[1]]
            type_pairs.append(temp)
        return type_pairs

    def count_contacts(self, contacts_types):

        unique_contacts = []
        count_array = []
        for con in contacts_types:
            if con in unique_contacts:
                count_array[unique_contacts.index(con)] += 1
            else:
                unique_contacts.append(con)
                count_array.append(1)
        return unique_contacts, count_array

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

    def read_itp_amino(self, itp_name):

        f = open(itp_name)
        lines = f.readlines()

        start = False
        tipes = []
        for line in lines:
            s = line.split()
            if len(s) == 0:
                start = False
            elif s[0].isnumeric() and start:
                tipes.append(s[3])
            elif s[:3] == ["[", "atoms", "]"]:
                start = True
        return tipes


    def contact_array(self, no_doubles=False, frame=0):

        time = self.trajectory[frame].time
        contacts = [[] for i in range(len(self.chains))]

        for i in range(len(self.chains)):
            for j in range(no_doubles * i, len(self.chains)):
                if self.check_contact(i, j) and i != j:
                    contacts[i].append(j)
        return contacts

    def combine_atom_groups(self, list_o_groups):

        indexes = []
        for group in list_o_groups:
            indexes.extend(group.indices)

        return self.universe.atoms[indexes]

    def make_graph(self, frame=0):

        G = nx.Graph()
        color_map = []

        for i in range(len(self.chains)):
            G.add_node(i)
            color_map.append('blue')

        contact_array = self.contact_array(no_doubles=True, frame=frame)
        for ind, line in enumerate(contact_array):
            for thing in line:
                G.add_edge(ind, thing)
        return G

    def get_clusters(self, frame=0):

        G = self.make_graph(frame=frame)

        S = [list(G.subgraph(c).copy().nodes) for c in nx.connected_components(G)]

        return S

    def get_all_cluster_shapes(self, frame=0, protein=True):

        clusters = self.get_clusters(frame=frame)
        aspheres = []
        for cluster in clusters:
            group = self.combine_atom_groups([self.chains[c] for c in cluster
                                              if len(np.unique(self.map_chains[c].resnames)) < 2 + 50 * protein])
            aspheres.append([self.group_asphericity(group, frame=frame),
                             self.group_acylindricity(group, frame=frame)])
        return aspheres


    def average_protein_chains_contacts_groups(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types = []
        polar_count = 0
        charged_count = 0
        non_polar_count = 0
        for frame in frames:
            print(frame)
            for protein_ind, protein_chain in enumerate(protein_chain_indices):
                if protein_chain in protein_chain_indices:
                    for chain_ind in group_chain_indices:
                        con = self.get_all_contacts(protein_chain, chain_ind, frame=frame)
                        types.extend(self.contacts_types(con, protein_chain, chain_ind, itp1=protein_itps[protein_ind]))
        typed, counted = self.count_contacts(types)
        # for i in range(len(typed)):
        # print(typed[i], counted[i])
        for ind in range(len(typed)):
            if "P" in typed[ind][0]:
                polar_count += counted[ind]
            elif "Q" in typed[ind][0]:
                charged_count += counted[ind]
            else:
                non_polar_count += counted[ind]

        return non_polar_count / len(frames), polar_count / len(frames), charged_count / len(frames)


    def average_protein_chains_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types = []

        #print(len(self.chains))
        #print(self.chains[0])
        #print(self.chains[1])
        #print(self.chains[2])
        #print("protein_chain_indices", protein_chain_indices)
        #print("group_chain_indices", group_chain_indices)
        #print("protein_itps", protein_itps)
        for frame in frames:
            print(frame)
            for protein_ind, protein_chain in enumerate(protein_chain_indices):
                if protein_chain in protein_chain_indices:
                    for chain_ind in group_chain_indices:
                        #print(protein_chain, chain_ind)
                        con = self.get_all_contacts(protein_chain, chain_ind, frame=frame)
                        types.extend(self.contacts_types(con, protein_chain, chain_ind, itp1=protein_itps[protein_ind]))
        typed, counted = self.count_contacts(types)

        types_zero = sorted(list(np.unique([typed[thing_ind][0] for thing_ind in range(len(typed))])))
        types_one = sorted(list(np.unique([typed[thing_ind][1] for thing_ind in range(len(typed))])))

        matrix = np.zeros((len(types_zero), len(types_one)))


        for ind in range(len(typed)):
            ind_zero = types_zero.index(typed[ind][0])
            ind_one = types_one.index(typed[ind][1])
            matrix[ind_zero][ind_one] = counted[ind]

        matrix = np.multiply(matrix, 1/(len(frames)))

        return types_zero, types_one, matrix


    def average_polymer_monomer_contacts(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types = []
        contacts_by_monomer = []
        self.init_mma_monomer_lists()


        for chain_ind in group_chain_indices:
            chain_contacts = np.zeros(len(self.chains_monomer_list[chain_ind]))
            #print(protein_chain, chain_ind)
            for protein_ind, protein_chain in enumerate(protein_chain_indices):
                if protein_chain in protein_chain_indices:
                    for frame in frames:
                        con = self.get_all_contacts(protein_chain, chain_ind, frame=frame)
                        #print(con)
                        monomers = self.contacts_monomer(con, protein_chain, chain_ind, itp1=protein_itps[protein_ind])
                        #print(monomers)
                        for thing in monomers:
                            chain_contacts[thing[1][0]] += 1
            contacts_by_monomer.append(chain_contacts)

        matrix = np.multiply(contacts_by_monomer, 1/(len(frames)))

        return matrix

    def reduced_average_protein_chains_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types_zero, types_one, matrix = self.average_protein_chains_contacts_matrix(frames, protein_chain_indices, group_chain_indices, protein_itps)

        new_matrix = np.zeros((2,2))

        new_types_zero = ["Hydrophobic", "Polar"]

        new_types_one = new_types_zero

        print(matrix)

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                print(types_zero[i], self.is_polar(types_zero[i]))
                print(types_one[j], self.is_polar(types_one[j]))
                new_matrix[self.is_polar(types_zero[i])][self.is_polar(types_one[j])] += matrix[i][j]
        return new_types_zero, new_types_one, new_matrix

    def is_polar(self, type):

        if type[0] == "C" or type[0] == "N":
            return 0
        if type[0] == "Q" or type[0] == "P":
            return 1
        #if type[0] == "C":
        #   return 0
        #if type[0] == "Q" or type[0] == "P" or type[0] == "N":
        #    return 1
        if type == "EO":
            return 1
        if type == "STY":
            return 0

        if type[0] == "S":
            return self.is_polar(type[1:])


    def is_polar_amino(self, AA):

        non_polar = ["ALA", "LEU", "VAL", "ILE", "PRO", "PHE", "MET", "TRP"]
        return AA not in non_polar


    def display_average_protein_chains_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types_zero, types_one, matrix = self.average_protein_chains_contacts_matrix(frames, protein_chain_indices, group_chain_indices, protein_itps)
        fig, ax = plt.subplots()
        im = ax.imshow(matrix, vmin=0, vmax=100)

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(types_one)))
        ax.set_yticks(np.arange(len(types_zero)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(types_one)
        ax.set_yticklabels(types_zero)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

        #Loop over data dimensions and create text annotations.
        for i in range(len(types_zero)):
            for j in range(len(types_one)):
                text = ax.text(j, i, str(int(matrix[i, j])),
                       ha="center", va="center", color="w", fontsize=10)

        print("sum of the matix: " + str(np.sum(matrix)))
        ax.set_title("Contact Matrix")
        fig.tight_layout()
        plt.savefig("contact_matrix.png")
        plt.show()

    def display_reduced_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):

        types_zero, types_one, matrix = self.reduced_average_protein_chains_contacts_matrix(frames, protein_chain_indices, group_chain_indices, protein_itps)
        fig, ax = plt.subplots()
        im = ax.imshow(matrix, vmin=0, vmax=100)

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(types_one)))
        ax.set_yticks(np.arange(len(types_zero)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(types_one)
        ax.set_yticklabels(types_zero)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

        #Loop over data dimensions and create text annotations.
        for i in range(len(types_zero)):
            for j in range(len(types_one)):
                text = ax.text(j, i, str(int(matrix[i, j])),
                       ha="center", va="center", color="w", fontsize=10)

        print("sum of the matix: " + str(np.sum(matrix)))
        ax.set_title("Contact Matrix")
        fig.tight_layout()
        plt.savefig("reduced_contact_matrix.png")
        plt.show()


    def average_protein_contact_groups(self, frames, protein_index, protein_itp):

        types = []
        polar_count = 0
        charged_count = 0
        non_polar_count = 0
        for frame in frames:
            print(frame)
            for ind, chain in enumerate(self.chains):
                if ind != protein_index and ind >= self.num_protein_chains:
                    #print(self.num_protein_chains, ind, protein_index)
                    con = self.get_all_contacts(protein_index, ind, frame=frame)
                    types.extend(self.contacts_types(con, protein_index, ind, itp1=protein_itp))
        typed, counted = self.count_contacts(types)
        #for i in range(len(typed)):
            #print(typed[i], counted[i])
        for ind in range(len(typed)):
            if "P" in typed[ind][0]:
                polar_count += counted[ind]
            elif "Q" in typed[ind][0]:
                charged_count += counted[ind]
            else:
                non_polar_count += counted[ind]

        return non_polar_count/len(frames), polar_count/len(frames), charged_count/len(frames)

    def graph_protein_contacts_time(self, frames, protein_index, protein_itp, save_name="contacts_time",
                                    time_unit =1, time_label="Time", chain_indices=None):

        non_polars = []
        polars = []
        charges = []
        totals = []
        if not type(protein_index) == list:
            protein_index = [protein_index]
            protein_itp = [protein_itp]
        for frame in frames:
            if chain_indices is None:
                chain_indices = [chain for ind, chain in enumerate(self.chains) if ind not in protein_index]
            #print(protein_index, chain_indices)
            #quit()
            n, p, q = self.average_protein_chains_contacts_groups([frame], protein_index, chain_indices, protein_itp)
            t = n + p + q
            non_polars.append(n)
            polars.append(p)
            charges.append(q)
            totals.append(t)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Protein Contacts by Type")
        ax1.set_xlabel(time_label)
        ax1.set_ylabel("Number of Contacts")
        ax1.plot(np.multiply(frames, time_unit), non_polars, label="Nonpolar")
        ax1.plot(np.multiply(frames, time_unit), polars, label="Polar")
        ax1.plot(np.multiply(frames, time_unit), charges, label="Charged")
        ax1.plot(np.multiply(frames, time_unit), totals, label="Total")

        plt.legend(fontsize=12)
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def center_frame(self, group, frame):

        time = self.trajectory[frame].time

        center = group.center_of_mass(pbc=True)
        dim = time.triclinic_dimensions
        box_center = np.sum(dim, axis=0) / 2
        self.trajectory.atoms.translate(box_center - center)

    def largest_cluster_eigs(self, frame=0):

        time = self.trajectory[frame]
        matrix = time.triclinic_dimensions

        directions = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]

        shifts = [np.matmul(matrix, directions[i]) for i in range(len(directions))]

        indices = [ind for chain in self.chains for ind in chain]

        original = []
        for chain in self.chains:
            original.extend(chain.positions)

        positions = cp.deepcopy(original)

        for shift in shifts:
            positions.extend(np.add(shift, original))

        no_doubles = True
        contact_array = [[] for i in range(len(positions))]
        dist_arr = distances.distance_array(np.array(positions), np.array(positions))
        cut = 5.3
        for i in range(len(positions)):
            for j in range(no_doubles * i, len(positions)):
                if dist_arr[i][j] < cut:
                    contact_array[i].append(j)

        G = nx.Graph()

        for i in range(len(positions)):
            G.add_node(i)

        for ind, line in enumerate(contact_array):
            for thing in line:
                G.add_edge(ind, thing)

        clusters = [list(G.subgraph(c).copy().nodes) for c in nx.connected_components(G)]

        numbers = [nx.number_of_nodes(clusters[i]) for i in range(len(clusters))]

        biggest = clusters[numbers.index(np.max(numbers))]

        pos_of_biggest = [positions[biggest] for thing in biggest]

        gyr_ten = self.positions_gyration_tensor(biggest)

        eigs = np.linalg.eigvals(gyr_ten)
        eigs.sort()
        print("eigs", eigs)

    def read_xvg(self, filename="rmsf.xvg"):

        f = open(filename, "r")
        data = f.readlines()
        labels = ['time']
        do = ""
        first = True
        for line in data:
                done = False
                s = line.split()
                if s[0][0] == "#":
                        done = True
                        #print("nothing doing", line)
                elif s[0][0] == "@" and not done:
                        if s[1][0] == "s":
                                l_strings = s[3:]
                                string = ""
                                for thing in l_strings:
                                        string += thing + " "
                                labels.append(string)
                                #print("added a label", string)
                                done = True
                        else:
                                done = True
                                #print("nothing doing", line)
                elif first and not done:
                        first = False
                        done = True
                        all_data = [[] for _ in range(len(labels))]
                        for ind, label in enumerate(labels):
                            all_data[ind].append(float(s[ind]))
                elif not done:
                        for ind, label in enumerate(labels):
                            all_data[ind].append(float(s[ind]))
        return all_data, labels

    def read_temperature(self, filename="temperature.xvg"):

        f = open(filename, 'r')

        data = f.readlines()

        labels = ['time']
        do = ""

        first = True

        for line in data:
            done = False
            s = line.split()
            if s[0][0] == "#":
                done = True
                print("nothing doing", line)
            elif s[0][0] == "@" and not done:
                if s[1][0] == "s":
                    l_strings = s[3:]
                    string = ""
                    for thing in l_strings:
                        string += thing + " "
                    labels.append(string)
                    print("added a label", string)
                    done = True
                else:
                    done = True
                    print("nothing doing", line)
            elif first and not done:
                first = False
                done = True
                all_data = [[] for _ in range(len(labels))]
                for ind, label in enumerate(labels):
                    all_data[ind].append(float(s[ind]))
            elif not done:
                for ind, label in enumerate(labels):
                    all_data[ind].append(float(s[ind]))

        return all_data, labels

    def read_rmsd(self, filename="rdf.xvg"):

          f = open(filename, 'r')

          data = f.readlines()

          labels = ['time']
          do = ""

          first = True


          for line in data:
                  done = False
                  s = line.split()
                  if s[0][0] == "#":
                          done = True
                          #print("nothing doing", line)
                  elif s[0][0] == "@" and not done:
                          if s[1][0] == "s":
                                  l_strings = s[3:]
                                  string = ""
                                  for thing in l_strings:
                                          string += thing + " "
                                  labels.append(string)
                                  #print("added a label", string)
                                  done = True
                          else:
                                  done = True
                                  #print("nothing doing", line)
                  elif first and not done:
                                  first = False
                                  done = True
                                  all_data = [[] for _ in range(len(labels))]
                                  for ind, label in enumerate(labels):
                                          all_data[ind].append(float(s[ind]))
                  elif not done:
                                  for ind, label in enumerate(labels):
                                          all_data[ind].append(float(s[ind]))
          return all_data, labels


    def graph_rmsf(self, filename="rmsf.xvg"):

        all_data, labels = self.read_xvg(filename=filename)
        save_name = filename[:-4] + ".png"
        self.graph_rmsf_data(all_data, labels, save_name=save_name)

    def graph_rmsf_data(self, all_data, labels, save_name="rmsf.png", ylim=[0,.5]):


         fig = plt.figure()
         #print(all_data)
         ax1 = fig.add_subplot(111)
         ax1.set_title("Backbone Root Mean Square Fluctuation")
         ax1.set_xlabel("Residue Number")
         ax1.set_ylabel('RMSF')
         for i in range(1, len(labels)):
             ax1.plot(np.add(29, range(len(all_data[i]))), all_data[i], label=labels[i])
         ax1.set_ylim(ylim)
         active = [88, 238, 159, 160, 185, 206, 237]
         #active = np.subtract(active, 29)
         for act in active:
             ax1.plot([act for _ in range(20)], [i * .05 for i in range(-10, 10)], label=str(act))
         plt.legend(fontsize=12)
         plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
         plt.show()

    def graph_rmsd_data(self, all_data, labels, filename="rdf.xvg", time_label="Time", time_unit=1, save_name="rdf.png"):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Backbone Root Mean Square Displacement")
        ax1.set_xlabel(time_label)
        ax1.set_ylabel("Root Mean Square Displacement (nm)")
        ax1.plot(np.multiply(range(len(all_data[1])), time_unit), all_data[1], label=labels[1])
        plt.legend(fontsize=12)
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_rmsd(self, filename="rdf.xvg"):

        all_data, labels = self.read_rmsd(filename=filename)
        save_name = filename[:-4] + ".png"
        self.graph_rmsd_data(all_data, labels, save_name=save_name)
        print("Average RMSD: ", np.average(all_data[1]))
        print("Standard Deviation: ", np.std(all_data[1]))

    def graph_contact_energy(self, filename="temperature.xvg"):

        if not type(filename) == type(["l", "i", "s", "t"]):
            all_data, labels = self.read_temperature(filename=filename)
            save_name = filename[:-4] + ".png"
            self.graph_contact_energy_data(all_data, labels, save_name=save_name, last=True)
        else:
            data = []
            for name in filename:
                print(name)
                self.graph_contact_energy(name)
                all_data, labels = self.read_temperature(filename=name)
                if len(all_data) > 2:
                    for i in range(2, len(all_data)):
                        all_data[1] = np.add(all_data[1], all_data[i])
                save_name = "combo.png"
                data.append(all_data[1])
            data = np.array(data).T
            data = [np.sum(row) for row in data]
            data = [all_data[0], data]
            self.graph_contact_energy_data(data, labels, save_name=save_name, last=True)


    def graph_contact_energy_data(self, all_data, labels,save_name="site_contact_energy.png", last=False):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Logged Energies")

        ax1.set_xlabel(labels[0])
        ax1.set_ylabel('KJ/mol')
        # all_data[2] = np.multiply(all_data[2],6/2)
        for i in range(1, len(all_data)):
            ax1.plot(all_data[0], all_data[i], label=labels[i])

            if last:
                print("Average Contact Energy: ", np.average(all_data[i][int(len(all_data[i]) *  3 / 4):]))
                print("Standard Deviation: ", np.std(all_data[i][int(len(all_data[i]) *  3 / 4):]),
                      "Standard Error: ", np.std(all_data[i][int(len(all_data[i]) *  3 / 4):]) / np.sqrt(len(all_data[i][int(len(all_data[i]) *  3 / 4):])))
            else:
                print("Average Contact Energy: ", np.average(all_data[i]))
                print("Standard Deviation: ", np.std(all_data[i]))

        plt.legend()
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def net_rmsf(self, original_file, new_file):

        all_data_o, labels_o = self.read_xvg(original_file)
        all_data_n, labels_n = self.read_xvg(new_file)

        all_data = np.subtract(all_data_n, all_data_o)


        print(np.sum(all_data[1][20:-5]), len(all_data[1][20:-5]), np.sum(all_data[1][20:-5]) / len(all_data[1][20:-5]) )

        self.graph_rmsf_data(all_data, labels_n, save_name="net_rmsf.png", ylim=[-.2, .2])

    def get_positions_chain(self, index, frame):

        time = self.trajectory[frame].time
        positions = []
        for ind, o_pos in enumerate(self.chains[index].positions):
            position = o_pos
            if ind != 0:
                position = self.new_pos(position, self.chains[index].positions[ind - 1])
            positions.append(position)

        return positions

    def new_pos(self, pos, ref, silence=False):

        do = False
        if la.norm(np.subtract(pos, ref)) > 50:
            #print("start")
            #print(pos, ref)
            do = True
        x_dist = pos[0] - ref[0]
        if x_dist > self.universe.dimensions[0] / 2:
            pos[0] = pos[0] - self.universe.dimensions[0]
        elif x_dist < - self.universe.dimensions[0]/2:
            pos[0] = pos[0] + self.universe.dimensions[0]

        y_dist = pos[1] - ref[1]
        if y_dist > self.universe.dimensions[1] / 2:
            pos[1] = pos[1] - self.universe.dimensions[1]
        elif y_dist < - self.universe.dimensions[1]/2:
            pos[1] = pos[1] + self.universe.dimensions[1]

        z_dist = pos[2] - ref[2]
        if z_dist > self.universe.dimensions[2] / 2:
            pos[2] = pos[2] - self.universe.dimensions[2]
        elif z_dist < - self.universe.dimensions[2]/2:
            pos[2] = pos[2] + self.universe.dimensions[2]
        if do:
         #   print(pos, ref)
         #   print(la.norm(np.subtract(pos, ref)))
            if la.norm(np.subtract(pos,ref)) > 50 and not silence:
                print("aaaaaaaa")
                print(pos, ref)
         #   print("done")
        return pos

    def water_touchers(self, frames, protein_index, protein_itp, cut=5.3):

        number = []
        others = self.universe.atoms[551:]
        petase = self.universe.atoms[:551]
        #itp = self.read_itp(protein_itp)
        itp = self.read_itp_amino(protein_itp)
        #print(itp)
        #quit()
        touchers = []
        nonpolars = []
        for frame in frames:
            print(frame)
            time = self.trajectory[frame].time
            dist_arr = distances.distance_array(petase.positions, others.positions,
                                            box=self.universe.dimensions)

            for i in range(dist_arr.shape[0]):
                for j in range(dist_arr.shape[1]):
                    if dist_arr[i][j] < cut:
                        touchers.append(i)
            number.append(len(np.unique(np.array(touchers))))

        touchers_total = []
        for i in range(551):
            if touchers.count(i) >= len(frames):
                touchers_total.append(i)

        print(touchers_total)
        nonpolars = np.sum([not self.is_polar_amino(itp[t]) for t in touchers_total])
        print(nonpolars/len(touchers_total))

        return len(touchers_total)

class JustGraphs(GromacsAnalysis):

    def __init__(self):
        print("just making graphs")

    def graph_pmf_windows(self, file_names):


        dist_sets = []
        for file_name in file_names:
            dist_sets.append(self.read_pmf_window(file_name))


        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Histograms")
        ax1.set_xlabel("Z (nm)")
        ax1.set_ylabel("Occurences")

        for ind, dists in enumerate(dist_sets):
            np_hist, np_be = np.histogram(dists)

            ax1.plot(np_be[:-1], np_hist, label=file_names[ind])

        plt.legend(fontsize=12)
        plt.savefig("wham_hists.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def read_pmf_window(self, file_name):

        dists = []
        f = open(file_name)
        data = f.readlines()
        for line in data[18:]:
            s =line.split()
            if len(s) > 1:
                dists.append(float(s[1]))

        return dists

    def graph_pmf(self, filename="profile.xvg", z_theta=False, z_theta_hex_hex=False):

        f = open(filename)
        data = f.readlines()
        com = [float(line.split()[0]) for line in data[18:]]
        pmf = [float(line.split()[1]) for line in data[18:]]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Protein Contacts by Type")
        ax1.set_xlabel("COM")
        ax1.set_ylabel("PMF")

        if z_theta or z_theta_hex_hex:
            params = [0.03265391, - 1.06512216, 7.05113996]
            if z_theta_hex_hex:
                params = [1.95262544e-05, 7.49336163e-04, -3.30998987e-03, -2.25091662e+00, 8.66784446e+00]
            x_averages = np.multiply(10, com)
            order = len(params) - 1
            com = [np.sum([params[i] * thing ** (order - i) for i in range(order + 1)]) for thing in x_averages]
            ax1.set_xlabel("$\Theta_{B}$")


        ax1.plot(com, pmf)

        plt.legend(fontsize=12)
        plt.savefig("wham_pmf.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def read_pull_file(self, name):

        time = []
        z = []
        f = open(name)
        data = f.readlines()
        for line in data[18:]:
            s = line.split()
            if len(s) > 1:
                z.append(float(s[1]))
                time.append(float(s[0]))

        return time, z

    def graph_force_extension(self, xname, fname, smooth=10, name_label=None):


        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('X(nm)')
        ax1.set_ylabel('F(pN)')
        time, x = self.read_pull_file(xname)
        time2, f = self.read_pull_file(fname)
        smooth_x = [np.average(x[smooth * i: smooth*i + 1]) for i in range(int(len(x)/smooth) - 1)]
        smooth_f = [np.average(f[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        smooth_f = np.multiply(smooth_f, 4.114 / 2.479)
        smooth_f = smooth_f[:min(len(smooth_f), len(smooth_x))]
        smooth_x = smooth_x[:min(len(smooth_f), len(smooth_x))]
        ax1.plot(smooth_x, smooth_f)
        #ax1.plot(x, f)
        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig("xtension.png", bbox_inches='tight', pad_inches=.2)
        if name_label is not None:
            plt.savefig("xtension"+ str(name_label)+".png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def pull_work(self, xname, fname, smooth=10, cut=50):

        time, x = self.read_pull_file(xname)
        time2, f = self.read_pull_file(fname)
        smooth_x = [np.average(x[smooth * i: smooth * i + 1]) for i in range(int(len(x) / smooth) - 1)]
        smooth_f = [np.average(f[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        smooth_f = np.multiply(smooth_f, 4.114 / 2.479)
        smooth_f = smooth_f[:min(len(smooth_f), len(smooth_x))]
        smooth_x = smooth_x[:min(len(smooth_f), len(smooth_x))]
        i = 0
        while smooth_x[i] < cut and len(smooth_x) > i+1:
            i += 1
        smooth_f = np.array(smooth_f[:i])
        smooth_x = np.array(smooth_x[:i])
        return np.trapz(smooth_f, smooth_x)

    def rupture_force(self, xname, fname, smooth=10, cut=50):

        time, x = self.read_pull_file(xname)
        time2, f = self.read_pull_file(fname)
        smooth_x = [np.average(x[smooth * i: smooth * i + 1]) for i in range(int(len(x) / smooth) - 1)]
        smooth_f = [np.average(f[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        smooth_f = np.multiply(smooth_f, 4.114 / 2.479)
        smooth_f = smooth_f[:min(len(smooth_f), len(smooth_x))]
        smooth_x = smooth_x[:min(len(smooth_f), len(smooth_x))]
        i = 0
        while smooth_x[i] < cut and len(smooth_x) > i+1:
            i += 1
        smooth_f = np.array(smooth_f[:i])

        return np.max(smooth_f)


    def graph_force_time(self, fname, smooth=10, name_label=None):


        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('time')
        ax1.set_ylabel('F(pN)')
        time2, f = self.read_pull_file(fname)
        smooth_f = [np.average(f[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        smooth_f = np.multiply(smooth_f, 4.114/2.479)
        smooth_time = [np.average(time2[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        ax1.plot(smooth_time, smooth_f)
        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig("forcetime.png", bbox_inches='tight', pad_inches=.2)
        if name_label is not None:
            plt.savefig("forcetime"+ str(name_label)+".png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_extension_time(self, xname, smooth=10, name_label=None):


        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('time')
        ax1.set_ylabel('X (nm)')
        time2, f = self.read_pull_file(xname)
        smooth_f = [np.average(f[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        smooth_time = [np.average(time2[smooth * i: smooth * i + 1]) for i in range(int(len(f) / smooth) - 1)]
        ax1.plot(smooth_time, smooth_f)
        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig("xtensiontime.png", bbox_inches='tight', pad_inches=.2)
        if name_label is not None:
            plt.savefig("xtensiontime"+ str(name_label)+".png", bbox_inches='tight', pad_inches=.2)
        plt.show()


    def calc_pmf(self, x_files, f_files, z_theta=False, write_name="pmf.txt", z_theta_hex_hex=False, z_theta_hex_tri=False,
                 z_theta_hex_tri_2=False, z_theta_tri_hex_2_del=False, z_theta_tri_hex_2_free=False, z_theta_j_j=False,
                 z_theta_a_j=False, z_theta_j_n=False, z_theta_j_n_low=False, z_theta_j_b1=False, z_theta_j_bdel=False,
                 z_theta_j_b2=False, z_theta_ho_hex_pent=False, z_theta_ho_hex_hex=False, z_theta_ho_hex_tri=False,
                 z_theta_ho_hex_tri2=False, z_theta_alpha_j=False, z_theta_alpha_alpha=False, z_theta_lower_ho=False,
                 z_theta_upper_ho=False, z_theta_eut=False, z_theta_csos=False, z_theta_EutMute=False,
                 auto_combine=False,
                 pent_sample_error=False, pent_pull_error=False, show_forces=False, pent_color=False, hex_error=False,
                 hex_color=False, hex_pull_error=False, b_color=False, b1_error=False, b2_error=False, pull_b_prime_error=False,
                 pull_b_free_error=False, alpha_alpha_pull_error=False, ho_hex_pent_error=False, ho_hex_hex_pull_error=False,
                 ho_hex_hex_error=False, jb1_error=False, jb1_pull_error=False, jb2_error=False, jb2_pull_error=False,
                 jbprime_error=False, jbprime_pull_error=False, ab1_error=False, ab2_error=False, abprime_error=False,
                 aj_error=False, aj_pull_error=False, jj_error=False, jj_pull_error=False, alpha_j_bend_error=False,
                alpha_j_pull_error=False, alpha_bend_error=False, j_erk_pull_error=False, j_n_low_pull_error=False,
                 j_n_flat_pull_error=False, ho_hex_pent_pull_error=False, ho_hex_tri_pull_error=False,
                 ho_hex_tri_error=False, ho_hex_tri2_error=False, ho_hex_tri2_pull_error=False,
                 ho_hex_hex_lower_bend_error=False, ho_hex_hex_upper_bend_error=False, ho_hex_hex_upper_pull_error=False,
                 eut_bend_error=False, csos_bend_error=False, eut_double_bend_error=False, eut_double_pull_error=False,
                 rmmh_pull_error=False, csos_pull_error=False, eut_pull_error=False):


        x_averages = []
        f_averages = []
        pmf = []
        if not auto_combine:
            for ind, name in enumerate(x_files):
                time, x = self.read_pull_file(name)
                #x = x[2*int(len(x) / 3):3*int(len(x) / 3)]
                x_averages.append(np.average(x))
                time, f = self.read_pull_file(f_files[ind])
                #f = f[2*int(len(f) / 3):3*int(len(f) / 3)]
                f_averages.append(np.average(f))
        else:
            x_averages = self.read_as_groups(x_files)
            f_averages = self.read_as_groups(f_files)

        #for i in range(len(x_averages)):
        #    print(x_averages[i], f_averages[i])

        temp = cp.deepcopy(x_averages)

        x_averages = sorted(x_averages, key=float)

        f_averages_2 = [f_averages[temp.index(thing)]  for thing in x_averages]

        f_averages = f_averages_2

        #print()
        #for i in range(len(x_averages)):
         #   print(x_averages[i], f_averages[i], f_averages_2[i])

        for ind in range(len(x_averages)):
            if ind == 0:
                pmf.append(0)
            else:
                slope = f_averages[ind-1]
                pmf.append(pmf[-1] + slope*(x_averages[ind] - x_averages[ind -1]))

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # ax1.set_title("Protein Contacts by Type")
        ax1.set_xlabel("Center of Mass Distance(nm)")
        if pent_color:
            ax1.set_xlabel("Y Distance, nm")
        ax1.set_ylabel("PMF (Kcal/mol)")

        units = 1 / 2.479
        units = units * .598

        if z_theta or z_theta_hex_hex or z_theta_hex_tri or z_theta_hex_tri_2 or z_theta_tri_hex_2_del or z_theta_tri_hex_2_free\
                or z_theta_j_j or z_theta_a_j or z_theta_j_n or  z_theta_j_n_low or z_theta_j_b1 or z_theta_j_bdel or z_theta_j_b2\
                or z_theta_ho_hex_pent or z_theta_ho_hex_hex or z_theta_ho_hex_tri or z_theta_ho_hex_tri2 or z_theta_alpha_j \
                or z_theta_alpha_alpha or z_theta_lower_ho or z_theta_upper_ho or z_theta_eut or z_theta_csos or z_theta_EutMute:
            params = [0.03265391, - 1.06512216, 7.05113996]
            if z_theta_hex_hex:
                params = [1.95262544e-05, 7.49336163e-04, -3.30998987e-03, -2.25091662e+00, 8.66784446e+00]
            if z_theta_alpha_alpha:
                params = [0.01000316, 1.6960565, -4.0500452]
            if z_theta_hex_tri:
                params = [-2.62643769e-05,  1.24201396e-03,  6.25559493e-02,  2.30684539e+00,  4.46811531e+01]
            if z_theta_hex_tri_2:
                params = [-1.11767460e-04, -3.69251928e-03, 3.37735408e-04,  2.57698832e+00, 5.04799186e+01]
            if z_theta_tri_hex_2_del:
                #params = [-2.72497177e-04, 2.02835210e-03, 6.08586898e-02, 1.69825177e+00, 2.47053902e+01]
                params = [-8.92650832e-05, - 9.35684515e-04, 3.62029754e-02, 2.01079259e+00, 2.51248418e+01]
            if z_theta_tri_hex_2_free:
                params = [ 1.96884677e-04,  3.99453187e-04, -1.18826842e-01,  1.93865813e+00, 6.12426686e+01]
            if z_theta_j_j:
                params = [ 0.03032936, -1.05645245,  0.20460532]
            if z_theta_a_j:
                #params = [0.01075255, -1.39031064, 5.8906275]
                params = [-0.01605193, - 1.92482842, 4.63249307]
            if z_theta_j_n:
                params = [0.05427288, -1.83987616, -1.87484278]
            if z_theta_j_n_low:
                params = [3.63840000e-02, 4.93699656e-01, 4.62156602e+01]
            if z_theta_j_b1:
                params = [1.68462150e-03, 2.47875558e+00, 1.58829498e+01]
            if z_theta_j_b2:
                #params = [-0.02301592, -1.07269287, 21.78050241]
                params = [-3.17241441e-04, - 1.26411509e+00, 2.10481780e+01]
            x_averages = np.multiply(-10, x_averages)
            if z_theta_j_bdel:
                params = [-0.03397273, 1.80729138, 12.94229078]
            if z_theta_alpha_j:
                params = [0.00306178, -1.83768295, -2.66792469]
            if z_theta_ho_hex_pent:
                params = [0.03327459, -1.1955925, -2.06676349]
            if z_theta_ho_hex_hex:
                params = [0.02765246, 1.16235195, -3.5457027]
            if z_theta_ho_hex_tri:
                params = [-0.0084523, - 1.57703713, 4.1062643]
            if z_theta_ho_hex_tri2:
                params = [0.01882712, -1.70688774, -1.56487244]
            if z_theta_lower_ho:
                params = [0.00736712, 1.7143699, - 2.55080351]
            if z_theta_upper_ho:
                params = [-5.58923122e-03, -2.15332733e+00, -5.61430473e+00]
            if z_theta_eut:
                params = [0.02535733, -1.15998358, 1.8297918]
            if z_theta_csos:
                params = [0.0080437, -1.63507024, 1.89130959]
            if z_theta_EutMute:
                params = [ 0.02976658, -1.46188215, -4.97187647]
            if z_theta or z_theta_hex_hex or z_theta_j_j or z_theta_a_j or z_theta_j_n or z_theta_j_n_low or z_theta_j_b1 \
                    or z_theta_j_bdel or z_theta_j_b2 or z_theta_ho_hex_pent or z_theta_ho_hex_tri or z_theta_ho_hex_tri2\
                    or z_theta_alpha_j or z_theta_alpha_alpha or z_theta_upper_ho or z_theta_eut or z_theta_csos or z_theta_EutMute:
                x_averages = np.multiply(-1, x_averages)
            #x_averages = [params[2] + params[1] * thing + params[0] * thing ** 2 for thing in x_averages]
            order = len(params) - 1
            x_averages = [np.sum([params[i] * thing ** (order - i) for i in range(order + 1)]) for thing in x_averages]
            ax1.set_xlabel(r"Bending Angle, $\theta_{B}$, ($^\circ$)")

        pent_sample_err = [2 - int(x > 25)*1.5 -int(np.abs(x-41.6) < .1)*.5 + .5 * int(x >60) for x in x_averages]

        pent_error = [np.abs((x-5.5)/1.5) + (.5* (int(x>5.35) + int(x>5.5))) for x in x_averages]
        hex_airer = [.5 + (int(x>0) * int(x < 34) * -.3 * x/34) + (int(x < 50) * int(x > 34) * np.abs(x-50)/16 * -.3)
                     + (1  * int(x> 50) * (x-50)/50) for x in x_averages]
        b_pull_1_error = [  2 * min(.5, x -6.3)   for x in x_averages]
        b_pull_2_error = [ min(.5, x -6.5) + int(x > 7) * 2 * min(.5, x-7)  for x in x_averages]

        pull_b_prime_airer = [ x- 6.5 for x in x_averages]

        pull_b_free_airer = [x - 6.5 * (2/3) for x in x_averages]

        hex_pull_airer = [ (x - 6.6) / .9 *2 for x in x_averages]

        ho_hex_pent_airer = [1 for x in x_averages]

        ho_hex_hex_airer = [1 for x in x_averages]

        jb1_airer = [ .5 + .5 * int(x>20) for x in x_averages]
        jb1_pull_airer = [ 5/4 * (x - min(x_averages)) for x in x_averages]

        jb2_airer = [max(0, x/30) for x in x_averages]
        jb2_pull_airer = [5/4 * (x - min(x_averages))  for x in x_averages]

        ab1_airer= [ np.abs(x-20) / 40 for x in x_averages]
        ab2_airer = [np.abs(x - 20) / 80 for x in x_averages]
        abprime_airer = [np.abs(x - 18) / 30 for x in x_averages]

        jprime_airer = [1 for x in x_averages]
        jprime_pull_airer = [5/4 * (x - min(x_averages)) for x in x_averages]

        alpha_alpha_pull_airer = [(x - 6.4)  for x in x_averages]
        ho_hex_hex_pull_airer = [2*(x - 6.4)  for x in x_averages]



        aj_airer= [ np.abs(x-18) * (int(x>18) /50 + int(x<18)/20) for x in x_averages]
        aj_pull_airer = [(x - min(x_averages)) * 1.5 for x in x_averages]

        jj_airer = [min(1, np.abs(x-25)/15) for x in x_averages]
        jj_pull_airer = [1.5 * (x - min(x_averages)) for x in x_averages]

        alpha_j_bend_airer = [min(1, np.abs(x-20)/30) for x in x_averages]

        alpha_bend_airer = alpha_j_bend_airer

        small_pull_airer = [.5 * (x - min(x_averages)) for x in x_averages]

        big_pull_airer = [2 * (x - min(x_averages)) for x in x_averages]

        ho_hex_tri_airer = [min(.5, np.abs(x - 30) / 15) for x in x_averages]

        ho_hex_tri2_airer = [.5 + max(0, x) / 20 for x in x_averages]

        ho_hex_hex_lower_bend_airer = [abs(x-10) * .5/10 for x in x_averages]
        ho_hex_hex_upper_bend_airer = [.2 + abs(x - 10) *  (int(x<10) *.1 / 10 + int(x>20) * 1.5/20) for x in x_averages]

        use = None
        if pent_sample_error:
            use = pent_sample_err
        if pent_pull_error:
            use = pent_error
        if hex_pull_error:
            use = hex_pull_airer
        if hex_error:
            use = hex_airer
        if b1_error:
            use = b_pull_1_error
        if b2_error:
            use = b_pull_2_error
        if pull_b_prime_error:
            use = pull_b_prime_airer
        if pull_b_free_error:
            use=pull_b_prime_airer
        if alpha_alpha_pull_error:
            use=alpha_alpha_pull_airer
        if ho_hex_pent_error:
            use = ho_hex_pent_airer
        if ho_hex_hex_pull_error:
            use = ho_hex_hex_pull_airer
        if ho_hex_hex_error:
            use = ho_hex_hex_airer
        if jb1_error:
            use=jb1_airer
        if jb1_pull_error:
            use = jb1_pull_airer
        if jb2_error:
            use = jb2_airer
        if jb2_pull_error:
            use = jb2_pull_airer
        if jbprime_error:
            use = jprime_airer
        if jbprime_pull_error:
            use = jprime_pull_airer
        if ab1_error:
            use = ab1_airer
        if ab2_error:
            use = ab2_airer
        if abprime_error:
            use = abprime_airer
        if aj_error:
            use=aj_airer
        if aj_pull_error:
            use = aj_pull_airer
        if jj_error:
            use = jj_airer
        if jj_pull_error:
            use = jj_pull_airer
        if alpha_j_bend_error:
            use = alpha_j_bend_airer
        if alpha_j_pull_error:
            use = jb1_pull_airer
        if alpha_bend_error:
            use = alpha_bend_airer
        if j_erk_pull_error:
            use = small_pull_airer
        if j_n_low_pull_error:
            use = small_pull_airer
        if j_n_flat_pull_error:
            use = jprime_pull_airer
        if ho_hex_pent_pull_error:
            use = jprime_pull_airer
        if ho_hex_tri_pull_error:
            use = jprime_pull_airer
        if ho_hex_tri_error:
            use=ho_hex_tri_airer
        if ho_hex_tri2_error:
            use = ho_hex_tri2_airer
        if ho_hex_tri2_pull_error:
            use = big_pull_airer
        if ho_hex_hex_lower_bend_error:
            use = ho_hex_hex_lower_bend_airer
        if ho_hex_hex_upper_bend_error:
            use = ho_hex_hex_upper_bend_airer
        if ho_hex_hex_upper_pull_error:
            use = [1.5 * (x - min(x_averages)) for x in x_averages]
        if eut_bend_error:
            use = [1.5 / 40 * np.abs(x - 40) * int(x<40) + .5/20 * np.abs(x-40) * int(x>40) for x in x_averages]
        if csos_bend_error:
            use = [1 / 25 * np.abs(x - 25) * int(x<25) + 1.5/45 * np.abs(x-25) * int(x>25) for x in x_averages]
        if eut_double_bend_error:
            use = [2 / 30 * np.abs(x - 30) * int(x<40) + 2/15 * np.abs(x-30) * int(x>30) for x in x_averages]
        if eut_double_pull_error:
            use = [2 * np.abs(x - min(x_averages)) for x in  x_averages]
        if rmmh_pull_error:
            use = [2 / .5 * np.abs(min(x - min(x_averages), .5)) for x in  x_averages]
        if csos_pull_error:
            use = [x - min(x_averages) for x in x_averages]
        if eut_pull_error:
            use = [x - min(x_averages) for x in x_averages]

        if not show_forces:
            if hex_color:
                ax1.errorbar(x_averages, np.multiply(pmf, units), yerr=use, color="#b1d2f2")
            elif b_color:
                ax1.errorbar(x_averages, np.multiply(pmf, units), yerr=use, color="red")
            #elif not pent_color:
             #   ax1.errorbar(x_averages, np.multiply(pmf, units), yerr=use)
            else:
                ax1.errorbar(x_averages, np.multiply(pmf, units), yerr=use, color="#6197d3")
        else:
            ax1.scatter(x_averages, f_averages)

        if write_name is not None:
            f = open(write_name, "w")
            print(use)
            print(pmf)
            for ind in range(len(x_averages)):
                f.write(str(x_averages[ind]) + " " + str(np.multiply(pmf, units)[ind]) + " " +
                        str(np.add(np.multiply(pmf, units)[ind], use[ind])) + " " + str(np.subtract(np.multiply(pmf, units)[ind], use[ind])) + "\n")
        #ax2 = fig.add_subplot(222)
        #ax2.plot(x_averages, f_averages)

        plt.legend(fontsize=12)
        #ax1.set_ylim((-1.1304961033631158, 12.884237570939225))
        print(x_averages, f_averages)
        plt.savefig("wham_pmf.png", bbox_inches='tight', pad_inches=.2)
        if pent_color or b_color:
            plt.savefig("wham_pmf.pdf", bbox_inches='tight', pad_inches=.2)
        plt.show()
        #print(ax1.get_ylim())


    def read_as_groups(self, file_names):

        groups = [[file_names[0]]]

        for name in file_names[1:]:
            found = False
            name_number = name.split("window")[1].split("_pull")[0]
            #print(name_number, name)
            for ind, group in enumerate(groups):
                group_number = group[0].split("window")[1].split("_pull")[0]
                #print(name_number, group_number, "group")
                word = "pull_down"
                #word = "pull_up"
                #word="extend"
                if (word in group[0] and word  in name) or (word not in group[0] and word not in name):
                    if group_number == name_number:
                        groups[ind].append(name)
                       #print("yeah", groups[ind])
                        found = True
            if not found:
                    groups.append([name])

        #print(groups)
        #quit()
        averages = []
        for group in groups:
            x_all = []
            for name in group:
                time, x = self.read_pull_file(name)
                x_all.extend(x)
            #x_all = x_all[1 * int(len(x) / 2):2 * int(len(x) / 2)]
            averages.append(np.average(x_all))
        return averages


    def dhs_model(self, velocity, k_0, k_m, x, k_s=1660):

        # 1660. pN/nm
        k_s = 208/41.14

        print(velocity, k_0, k_m, x, k_s)
        tau = x / velocity

        times = [tau * i /.01 for i in range(100)]
        dt = tau/100
        survive = [self.dhs_integral(velocity, k_s, k_m, x, k_0, t) for t in times]
        #quit()
        grated = np.trapz(survive, dx=dt)
        print(velocity, grated, velocity * grated, x)

        return -.04114 * k_s * (x - velocity * grated)

    def dhs_integral(self, v, k_s, k_m, x, k_0, t):

        top = k_0 * (np.exp(k_s * (-(x * x)/2) + v*x*t - .5 * (k_s * v * t) **2 / (k_m + k_s))- np.exp(-k_s*x*x/2) )
        bottom = v * k_s * x * np.power(k_m/ (k_m + k_s), 3/2)
        #print("time", t, "top bottom", top, bottom)
        return np.exp(-top /bottom)

    def simple_model(self, velocity, x, k):

        k_s = 1000 / 2.479 * 4.114
        return np.log(k_s * velocity * x  *np.exp(-.5772) * 1/k) / x


    def get_dhs_params(self, velocities, forces):

        guess= [1,1,1]
        parameters, covariance = curve_fit(self.dhs_model, velocities, forces, p0=guess)
        print("params", parameters)
        return parameters

    def get_simple_params(self, velocities, forces):

        guess = [1,1]
        parameters, covariance = curve_fit(self.simple_model, velocities, forces, p0=guess)
        print("params", parameters)
        return parameters

    def graph_dhs_model(self, velocities, forces):

        #forces = np.divide(forces, 4.114)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('Pull Velocity(A/ps)')
        ax1.set_ylabel('Average Rupture Force(pN)')
        #fit_k_0,fit_k_m, fit_x = self.get_dhs_params(velocities, forces)
        #fit_par, fit_par2 = self.get_simple_params(velocities, forces)
        #fit_par, fit_par2 = (.01,.0001)
        fit_vels = [np.power(10, i) for i in range(2, 13, 1)]
        fit_forces = [self.dhs_model(fit_vels[i], 1.6e-11, 286, 3.85, k_s=208.4) for i in range(len(fit_vels))]
        #fit_forces = [self.simple_model(fit_vels[i], fit_par, fit_par2) for i in range(len(fit_vels))]
        ax1.plot(np.divide(np.log(velocities), np.log(10)), forces, 'o', label="data")
        ax1.plot(np.divide(np.log(fit_vels), np.log(10)), fit_forces, '-', label="fit")
        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig("fit.png", bbox_inches='tight', pad_inches=.2)
        plt.show()



class GromacsAnalysisSet(object):

    def __init__(self, gro_name, xtc_name, directories, count=1, itps=None):

        self.directories = directories
        self.anals = [GromacsAnalysis(direct + "/" + gro_name, direct + "/" + xtc_name, count=count, itps=itps) for direct in directories]


    def graph_hist_z_one_polymer(self, chain_index, save_name="z_hist.png"):

        distances = []

        list_frames = [list(range(int(3 * len(a.trajectory) / 4), len(a.trajectory), 100))  for a in self.anals ]

        for ind, a in enumerate(self.anals):
            for frame in list_frames[ind]:
                distances.extend(list(a.z_dist_chain(chain_index, frame)))
                print(len(distances))


        bins = [4 * i for i in range(50)]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        #ax1.set_title("Frame " + str(frames[-1]))
        ax1.set_xlabel('z (nm)')
        ax1.set_ylabel('P(z)')




        hist, bin_edges = np.histogram(distances, bins=bins)
        bin_edges = np.divide(bin_edges, 10)
        print(hist, bin_edges)
        box_area = 400
        z_space = bin_edges[1] - bin_edges[0]
        z_space = 1
        num_frames = np.sum([len(frames) for frames in list_frames])
        total_factor = box_area *  z_space * num_frames
        hist = np.divide(hist, total_factor)
        #hist = np.divide(hist, np.sum(hist))
        ax1.plot(bin_edges[:-1], hist)

        np.savetxt("bin_edges.txt", bin_edges)
        np.savetxt("hist.txt", hist)

        plt.legend(fontsize=14, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()


    def count_protein_times_contacted(self, frames, cut=5.3, protein_index=0, chain_indices=None):

        count = np.zeros(len(self.anals[0].chains[protein_index].positions))


        for anal in self.anals:
            for frame in frames:
                count = np.add(count, anal.protein_times_contacted(protein_index, frame=frame, cut=cut, chain_indices=chain_indices))

        return np.divide(count, (len(frames)* len(self.anals)))



    def count_protein_contact_groups(self, frames, protein_index, protein_itp):


        if not type(protein_index) == list:
            protein_index = [protein_index]
            protein_itp = [protein_itp]

        polar_count = np.zeros((len(self.anals), len(frames)))
        charged_count = np.zeros((len(self.anals), len(frames)))
        non_polar_count = np.zeros((len(self.anals), len(frames)))
        for anal_ind, anal in enumerate(self.anals):
            for frame_ind, frame in enumerate(frames):
                print(frame)
                types = []
                #print(anal.num_protein_chains)
                #quit()
                for ind in list(range(anal.num_protein_chains, len(anal.map_chains))):
                        for protein_ind, protein in enumerate(protein_index):
                            con = anal.get_all_contacts(protein, ind, frame=frame)
                            types.extend(anal.contacts_types(con, protein, ind, itp1=protein_itp[protein_ind]))
                typed, counted = anal.count_contacts(types)

                for ind in range(len(typed)):
                    if "P" in typed[ind][0]:
                        polar_count[anal_ind][frame_ind] += counted[ind]
                    elif "Q" in typed[ind][0]:
                        charged_count[anal_ind][frame_ind] += counted[ind]
                    else:
                        non_polar_count[anal_ind][frame_ind] += counted[ind]

        return non_polar_count, polar_count, charged_count

    def graph_average_rmsf_data(self, filename="rmsf.xvg"):

        save_name = filename[:-4] + ".png"
        data, all_labels = self.anals[0].read_xvg(self.directories[0] + "/" + filename)
        all_data = [data[1]]
        for ind, anal in enumerate(self.anals[1:]):
            d, a_l = anal.read_xvg(self.directories[1 + ind] + "/" + filename)
            all_data.append(data[1])

        average_data = [[], np.average(all_data, axis=0)]
        self.anals[0].graph_rmsf_data(average_data, all_labels, save_name=save_name)


    def graph_all_rmsd_data(self, filename="rdf.xvg", chop=False):

        save_name = filename[:-4] + ".png"
        data, all_labels = self.anals[0].read_rmsd(self.directories[0] + "/" + filename)
        for ind, anal in enumerate(self.anals[1:]):
            d, a_l = anal.read_rmsd(self.directories[1+ind] + "/" + filename)
            if chop:
                d[0] = d[0][int(len(d[0]) / 2):]
                d[1] = d[1][int(len(d[1]) / 2):]
            data[0].extend(d[0])
            data[1].extend(d[1])
            #all_labels.extend(a_l)
        self.anals[0].graph_rmsd_data(data, all_labels, save_name=save_name)
        print("Average RMSD: ", np.average(data[1]))
        print("Standard Deviation: ", np.std(data[1]))
        return np.average(data[1]), np.std(data[1])

    def graph_all_contact_site_energy(self, filename="temperature.xvg"):

        save_name = filename[:-4] + ".png"
        data, all_labels = self.anals[0].read_rmsd(self.directories[0] + "/" + filename)
        for ind, anal in enumerate(self.anals[1:]):
            d, a_l = anal.read_rmsd(self.directories[1+ind] + "/" + filename)
            data[0].extend(np.add(d[0],data[0][-1]))
            data[1].extend(d[1])
            data[2].extend(d[2])
            data[3].extend(d[3])
            data[4].extend(d[4])

            #all_labels.extend(a_l)
        self.anals[0].graph_contact_energy_data(data, all_labels, save_name=save_name)
        print("Average Energy: ", np.average(np.sum(data[1:], axis=0)))
        print("Standard Deviation: ", np.std(np.sum(data[1:], axis=0)))

    def graph_all_contact_site_energy_dist(self, filename="temperature.xvg"):
        save_name = filename[:-4] + "_dist.png"
        data, all_labels = self.anals[0].read_rmsd(self.directories[0] + "/" + filename)
        for ind, anal in enumerate(self.anals[1:]):
            d, a_l = anal.read_rmsd(self.directories[1 + ind] + "/" + filename)
            data[0].extend(np.add(d[0], data[0][-1]))
            data[1].extend(d[1])
            data[2].extend(d[2])
            data[3].extend(d[3])
            data[4].extend(d[4])

        data = np.array(data)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Protein Active Site - Polymer Contact Energy")
        ax1.set_xlabel("Contact Energy (kJ/mol)")
        ax1.set_ylabel("P")

        d1_hist, d1_be = np.histogram(data[1].flatten())
        d2_hist, d2_be = np.histogram(data[2].flatten())
        d3_hist, d3_be = np.histogram(data[3].flatten())
        d4_hist, d4_be = np.histogram(data[4].flatten())

        ax1.plot(d1_be[1:], d1_hist, label=all_labels[1])
        ax1.plot(d2_be[1:], d2_hist, label=all_labels[2])
        ax1.plot(d3_be[1:], d3_hist, label=all_labels[3])
        ax1.plot(d4_be[1:], d4_hist, label=all_labels[4])

        plt.legend(fontsize=12)
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_protein_contacts_dist(self, frames, protein_index, protein_itp, save_name="contacts_dist"):


        n, p, q = self.count_protein_contact_groups(frames, protein_index, protein_itp)
        t = np.add(np.add(n, p), q)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Protein Contacts by Type")
        ax1.set_xlabel("Number of Contacts")
        ax1.set_ylabel("P")

        np_hist, np_be = np.histogram(n.flatten())
        p_hist, p_be = np.histogram(p.flatten())
        c_hist, c_be = np.histogram(q.flatten())
        t_hist, t_be = np.histogram(t.flatten())

        ax1.plot(np_be[:-1], np_hist, label="Nonpolar")
        ax1.plot(p_be[:-1], p_hist, label="Polar")
        ax1.plot(c_be[:-1], c_hist, label="Charged")
        ax1.plot(t_be[:-1], t_hist, label="Total")

        plt.legend(fontsize=12)
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_spatial_contacts(self, frames, protein_index=0, cut=5.3, save_name="contact_spatial.png", petase=True, p450=False, chain_indices=None):

        cvs_o = self.count_protein_times_contacted(frames, cut=cut, protein_index=protein_index, chain_indices=chain_indices)
        positions = self.anals[int(len(self.anals) - 1)].map_chains[protein_index].positions
        com = np.average(positions, axis=0)
        positions = np.subtract(positions, com)
        quat = QuaternionBetween(positions[118], [1, 0, 0])
        for ind, position in enumerate(positions):
            positions[ind] = quat.orient(position)

        w_spherical = self.anals[0].appendSpherical_np(positions)
        phi_o = w_spherical[:, 4]
        theta_o = w_spherical[:, 5]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_title("Spatial Contact Distribution")

        ax1.set_xlabel('Theta')
        ax1.set_ylabel('Phi')
        to_remove = []
        for ind in range(len(cvs_o)):
            if cvs_o[ind] < 1:
                to_remove.append(ind)
        phi = np.delete(phi_o, to_remove)
        theta = np.delete(theta_o, to_remove)
        cvs = np.delete(cvs_o, to_remove)
        print(len(phi))

        active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                        423, 424, 425, 426, 427, 428]

        f_g_loop = list(range(399, 456, 1))
        ax1.set_xlim(-3.14, 3.14)
        ax1.set_ylim(0, 3.14)
        ax1.scatter(theta, phi, c=cvs, cmap=cm.viridis, vmin=1, vmax=5)
        if petase:
            for act in active_site:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act], c="red")
                else:
                    ax1.scatter(theta_o[act], phi_o[act], c="red", marker="x")
        if p450:
            for act in f_g_loop:
                if act in to_remove:
                    ax1.scatter(theta_o[act], phi_o[act], c="red")
                else:
                    ax1.scatter(theta_o[act], phi_o[act], c="red", marker="x")

        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_contacts(self, frames, protein_index=0, cut=5.3, save_name="contact_numbers.png", petase=True, p450=False, write=True, chain_indices=None):

            cvs = self.count_protein_times_contacted(frames, cut=cut, protein_index=protein_index, chain_indices=chain_indices)

            fig = plt.figure()

            ax1 = fig.add_subplot(111)

            ax1.set_title("Average Number of Contacts")

            active_site = [117, 118, 267, 268, 269, 270, 271, 272, 312, 313, 314, 315, 316, 360, 361,
                           423, 424, 425, 426, 427, 428]


            f_g_loop = list(range(399,456,1))


            ax1.set_xlabel('Bead index')
            ax1.set_ylabel('Average number of contacts')
            ax1.plot(list(range(len(cvs))), cvs)
            if petase:
                for act in active_site:
                    ax1.plot([act - 1 for _ in range(10)], [i * .5 for i in range(10)],color="red", label="Active Site")
            if p450:
                for act in f_g_loop:
                    ax1.plot([act - 1 for _ in range(10)], [i * .5 for i in range(10)],color="red", label="F-G Channel")

            if write:
                f = open("protein_index_contact_numbers.txt", 'w')
                for i in range(len(cvs)):
                    f.write(str(i) + " " + str(cvs[i]) + "\n")
            #plt.legend()
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
            plt.show()

    def graph_petase_hist(self, index, frames):


        unique_names, rdf_distances = self.anals[0].rdf_polymers_petase_center(index, frames)
        for anal in self.anals[1:]:
            unique_names_1, rdf_distances_1 = anal.rdf_polymers_petase_center(index, frames)
            unique_names, rdf_distances = self.matrix_add_1d(unique_names, rdf_distances, unique_names_1, rdf_distances_1)


        bins = [1 * i for i in range(0, 180, 5)]
        fig = plt.figure()

        print(unique_names)
        #print(distances)


        #combined = cp.deepcopy(rdf_distances[0])
        #combined.extend(rdf_distances[2])
        #combined.extend(rdf_distances[4])

        #unique_names.append("combined")
        #rdf_distances.append(combined)

        #combined2 = cp.deepcopy(rdf_distances[1])
        #combined2.extend(rdf_distances[5])

        #unique_names.append("combined2")
        #rdf_distances.append(combined2)

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('R (nm)')
        ax1.set_ylabel('P(R)')

        myorder = [6, 7, 3, 5, 1, 0, 4, 2]

        np_list = [0, 2, 4, 6]

        real_names = ["C1", "EO", "Na", "Qa", "SC1", "SP2", "combined", "combined2"]
        names = ["NonPolar (C1)", "Polar (EO)", "NonPolar (Na)", "Negatively Charged (Qa)", "NonPolar (SC1)", "Polar (SP2)", "NonPolar (All)", "Polar (All)"]
        colors = ["#aaaa00",  "#0000ff", "#aaaa00", "#00ff7f", "#aaaa00",  "#0000ff", "#aaaa00", "#0000ff"]
        shapes = ["dashed", "dashed", "dotted", "solid", "dashdot", "dotted", "solid", "solid"]

        #for ind in myorder:
        for ind in range(len(unique_names)):
            name_ind = real_names.index(unique_names[ind])
            hist, bin_edges = np.histogram(rdf_distances[ind], bins=bins)
            #if ind in np_list:
                #hist = np.divide(hist, 2)
            hist = np.divide(hist, np.sum(hist))
            bin_edges = np.divide(bin_edges, 10)
            #if shapes[ind] == "solid":
             #   mpl.rcParams['lines.linewidth'] = 2
            #else:
            #    mpl.rcParams['lines.linewidth'] = 1
            ax1.plot(bin_edges[:-1], hist, label=names[name_ind], color=colors[name_ind], linestyle=shapes[name_ind])

        plt.legend(fontsize=10, loc="upper right", ncol=1, bbox_to_anchor=(1.1,1.0))
        plt.savefig("petase_hist.png", bbox_inches='tight', pad_inches=.2)
        plt.show()



    def display_average_protein_chains_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):


        types_zero, types_one, matrix = self.anals[0].average_protein_chains_contacts_matrix(frames, protein_chain_indices,
                                                                                    group_chain_indices, protein_itps)

        for anal in self.anals[1:]:
            types_zerop, types_onep, matrixp = anal.average_protein_chains_contacts_matrix(frames, protein_chain_indices, group_chain_indices, protein_itps)
            types_zero, types_one, matrix = self.matrix_add(types_zero, types_one, matrix, types_zerop, types_onep, matrixp)


        matrix = np.multiply(matrix, 1/ len(self.anals))
        fig, ax = plt.subplots()
        im = ax.imshow(matrix, vmin=0, vmax=100)

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(types_one)))
        ax.set_yticks(np.arange(len(types_zero)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(types_one)
        ax.set_yticklabels(types_zero)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="anchor", fontsize=14)

        #Loop over data dimensions and create text annotations.
        for i in range(len(types_zero)):
            for j in range(len(types_one)):
                text = ax.text(j, i, str(int(matrix[i, j])),
                       ha="center", va="center", color="w", fontsize=10)

        print("sum of the matix: " + str(np.sum(matrix)))
        ax.set_title("Contact Matrix")
        fig.tight_layout()
        plt.savefig("contact_matrix.png")
        plt.show()


    def display_reduced_contacts_matrix(self, frames, protein_chain_indices, group_chain_indices, protein_itps):


        types_zero, types_one, matrix = self.anals[0].reduced_average_protein_chains_contacts_matrix(frames, protein_chain_indices,
                                                                                    group_chain_indices, protein_itps)

        for anal in self.anals[1:]:
            types_zerop, types_onep, matrixp = anal.reduced_average_protein_chains_contacts_matrix(frames, protein_chain_indices, group_chain_indices, protein_itps)
            types_zero, types_one, matrix = self.matrix_add(types_zero, types_one, matrix, types_zerop, types_onep, matrixp)


        matrix = np.multiply(matrix, 1/ len(self.anals))
        fig, ax = plt.subplots()
        im = ax.imshow(matrix, vmin=0, vmax=100)

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(types_one)))
        ax.set_yticks(np.arange(len(types_zero)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(types_one)
        ax.set_yticklabels(types_zero)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="anchor", fontsize=14)

        #Loop over data dimensions and create text annotations.
        for i in range(len(types_zero)):
            for j in range(len(types_one)):
                text = ax.text(j, i, str(int(matrix[i, j])),
                       ha="center", va="center", color="w", fontsize=10)

        print("sum of the matix: " + str(np.sum(matrix)))
        ax.set_title("Contact Matrix")
        fig.tight_layout()
        plt.savefig("reduced_contact_matrix.png")
        plt.show()


    def display_polymer_contacts(self, frames, protein_chain_indices, group_chain_indices, protein_itps):


        matrix = self.anals[0].average_polymer_monomer_contacts(frames, protein_chain_indices,
                                                                                    group_chain_indices, protein_itps)

        types = self.anals[0].chains_monomer_list


        for anal in self.anals[1:]:
            new_matrix = anal.average_polymer_monomer_contacts(frames, protein_chain_indices,
                                                                    group_chain_indices, protein_itps)
            for chain_ind in range(len(new_matrix)):
                matrix[chain_ind] = matrix[chain_ind] + new_matrix[chain_ind]
            print(matrix)


        #quit()

        for ind, chain in enumerate(matrix):

            convert_types = []
            for t in types[ind + self.anals[0].num_protein_chains]:
                #print("t", t)
                if t == "EHMA":
                    convert_types.append("H")
                elif t == "SPMA":
                    convert_types.append("-")
                else:
                    convert_types.append(" ")
            #print(types)
            #print(types[ind + self.anals[0].num_protein_chains])
            #print(convert_types)
            #quit()
            chain = np.multiply(chain, 1 / len(self.anals))
            fig, ax = plt.subplots()
            ax.plot(chain)

            # We want to show all ticks...
            ax.set_xticks(np.arange(len(convert_types)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(convert_types)

            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor", fontsize=8)
            print(str(convert_types.count("H")) + " hydrophobic monomers in chain " + str(ind + 1))
            print(str(convert_types.count("-")) + " charged monomers in chain " + str(ind + 1))
            print("chain sum: " + str(np.sum(chain)))
            print("sum of the matix: " + str(np.sum(matrix)/len(self.anals)))
            ax.set_title("Contact Matrix")
            fig.tight_layout()
            self.write_monomer_contact_file(convert_types, chain, "chain_"+str(ind + 1)+".txt" )
            plt.savefig("chain_"+str(ind + 1)+"_matrix.png")
            plt.show()

    def write_monomer_contact_file(self, convert_types, chain, name):

        f = open(name, 'w')

        for i in range(len(convert_types)):
            if convert_types[i] == " ":
                convert_types[i] = "L"

        for i in range(len(chain)):
            f.write(convert_types[i] + " " + str(chain[i]) + "\n")
        f.close()

    def matrix_add(self, types_00, types_01, matrix_0, types_10, types_11, matrix_1):

        new_types0 = sorted(list(np.unique(types_00 + types_10)))
        new_types1 = sorted(list(np.unique(types_01 + types_11)))
        new_matrix = np.zeros((len(new_types0), len(new_types1)))

        for ind, name in enumerate(types_00):
            for ind2, name2 in enumerate(types_01):
                new_index0 = new_types0.index(name)
                new_index1 = new_types1.index(name2)
                new_matrix[new_index0][new_index1] = new_matrix[new_index0][new_index1] + matrix_0[ind][ind2]

        for ind, name in enumerate(types_10):
            for ind2, name2 in enumerate(types_11):
                new_index0 = new_types0.index(name)
                new_index1 = new_types1.index(name2)
                new_matrix[new_index0][new_index1] = new_matrix[new_index0][new_index1] + matrix_1[ind][ind2]

        return new_types0, new_types1, new_matrix

    def matrix_add_1d(self, names_0, data_0, names_1, data_1):

        names_new = cp.deepcopy(names_0)
        data_new = cp.deepcopy(data_0)

        for ind, name in enumerate(names_1):
            if name in names_new:
                put_index = names_new.index(name)
                data_new[put_index].extend(data_1[ind])
            else:
                data_new.append(data_1[ind])
                names_new.append(names_1[ind])

        return names_new, data_new


class JustGraphsSet(GromacsAnalysisSet):

    def __init__(self, directories):

        self.directories = directories
        self.anals = [JustGraphs() for direct in directories]