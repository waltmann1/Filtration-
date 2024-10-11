from __future__ import division
import gsd.hoomd
import gsd.fl
import numpy as np
import numpy.linalg as la
import os.path
import networkx as nx
import copy as cp
from matplotlib import pyplot as plt
from Analysis import Analysis
plt.rcParams.update({'font.size': 22})
import hoomd.data as hbox
from mpl_toolkits.mplot3d import Axes3D
import math as m
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import rms


class NetworkAnalysis(object):

    def __init__(self, gsd_name, map_name, protein_pdb=None, protein_itp=None):
        f = gsd.fl.open(name=gsd_name, mode='rb', application='', schema='hoomd',
                        schema_version=[1, 0])
        self.trajectory = gsd.hoomd.HOOMDTrajectory(f)
        self.dyes = None
        self.chain_indices = None
        self.chain_sequences = None
        self.all_chain_indices = None
        self.read_map(map_name)
        self.frames = []
        self.box = self.trajectory.read_frame(0).configuration.box[:3]
        self.gsd_name = gsd_name

        self.go_protein = None
        if protein_itp is not None and protein_pdb is not None:
            from GoProteinDye import GoProteinDye
            self.go_protein = GoProteinDye(protein_pdb, protein_itp)


    def get_box(self, frame):
        self.box = self.trajectory.read_frame(frame).configuration.box[:3]


    def read_map(self, map_name):

        f = open(map_name, 'r')

        dyes = []
        chain_indices = []
        chain_sequences = []

        data = f.readlines()

        for line in data:
            s = line.split()
            if s[0] == "dye":
                dyes.append([int(float(s[i])) for i in range(1, len(s))])
            elif s[0] == "rando_polymer":
                i = 2
                mon = []
                chain = []
                sequence = [s[1]]
                while i < len(s):
                    if not any(char.isalpha() for char in s[i]):
                        mon.append(int(float(s[i])))
                    else:
                        chain.append(mon)
                        mon = []
                        sequence.append(s[i])
                    i += 1
                chain.append(mon)
                chain_indices.append(chain)
                chain_sequences.append(sequence)

        self.dyes = dyes
        self.chain_sequences = chain_sequences
        self.chain_indices = chain_indices
        self.all_chain_indices = [[self.dyes[i]] for i in range(len(self.dyes))] + self.chain_indices

    def get_position(self, index, frame):

        return frame.particles.position[index]

    def get_position_min_image(self, index, frame):

        #print(frame.particles.image[index], frame.configuration.box[:3])
        return frame.particles.position[index] + np.multiply(frame.particles.image[index], frame.configuration.box[:3])

    def get_type(self, index, frame):

        if not isinstance(index, int):
            # print("Warning index "  + str(index) + "is not an integer")
            return ("!!!Faketype")
        return frame.particles.types[frame.particles.typeid[index]]


    def look_for_connection(self, chain1, chain2, frame, cut=0.5):

        #print("checking new")
        positions1 = np.array([self.get_position(index1, frame) for mon1 in chain1 for index1 in mon1
                      if self.get_type(index1, frame) == "A" or len(chain1) == 1])
        positions2 = np.array([self.get_position(index2, frame) for mon2 in chain2 for index2 in mon2
                      if self.get_type(index2, frame) == "A" or len(chain2) == 1])

        dist_array = distances.distance_array(positions1, positions2)
        if np.any(dist_array < cut):
            return True

    def look_for_connections_all(self, not_in_cluster,current_cluster,current_index, frame):

        current_cluster.append(current_index)
        not_in_cluster.remove(current_index)
        for chain_index in not_in_cluster:
            if self.look_for_connection(self.all_chain_indices[current_index], self.all_chain_indices[chain_index], frame):
                print("hit", current_index, chain_index)
                self.look_for_connections_all(not_in_cluster, current_cluster, chain_index, frame)

    def get_clusters_all(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters2.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters2.txt")

        not_in_cluster = list(range(0, len(self.all_chain_indices)))

        clusters = []

        current_cluster = []

        while len(not_in_cluster) != 0:
            print("not in cluster", not_in_cluster)
            self.look_for_connections_all(not_in_cluster, current_cluster, not_in_cluster[0], frame)
            clusters.append(current_cluster)
            current_cluster = []

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters2.txt", "w")
        for cluster in clusters:
            string = ""
            for thing in cluster:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return clusters


    def read_data(self, path, seq=False):

        f = open(path)
        data = f.readlines()
        out = []
        for line in data:
            s = line.split()
            if not seq:
                part = [int(thing) for thing in s]
            else:
                part = [str(thing) for thing in s]
            out.append(part)
        return out

    def get_all_connections_from_clusters(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters4.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters4.txt")


        clusters = []

        current_cluster = []

        known = self.get_clusters_all(number)

        for index1, chain_1 in enumerate(self.all_chain_indices):
            current_cluster.append(index1)
            for cluster_index, known_cluster in enumerate(known):
                if index1 in known_cluster:
                    for index_2 in known_cluster:
                        chain_2 = self.all_chain_indices[index_2]
                        if self.look_for_connection(chain_1, chain_2, frame) and index1 != index_2:
                            current_cluster.append(index_2)
                            print("current cluster", current_cluster)
                            print("hit " + str(index1) + " " + str(index_2))
            clusters.append(current_cluster)
            print("clusters", clusters)
            current_cluster = []
            #print("clusters ", clusters)

        #while len(not_in_cluster) != 0:
        #    print("not in cluster", not_in_cluster)
        #    self.look_for_connections_charged(self.all_chain_indices, frame)
        #    clusters.append(current_cluster)
         #   current_cluster = []
        print("writing")
        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters4.txt", "w")
        for cluster in clusters:
            print(cluster)
            string = ""
            for thing in cluster:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return clusters

    def make_graph(self, frame, include_dye=True):

        G = nx.Graph()
        color_map = []
        total_number = len(self.all_chain_indices)
        dye_number = len(self.dyes)
        chain_number = total_number - dye_number
        for i in range(total_number):
            G.add_node(i)
            if i < dye_number:
                color_map.append('blue')
            else:
                color_map.append('r')

        clusters = self.get_all_connections_from_clusters(frame)

        for i in range( total_number):
            for thing in clusters[i]:
                    if i != thing:
                        G.add_edge(i, thing)

        if not include_dye:
            G.remove_nodes_from(range(len(self.dyes)))
            color_map = color_map[len(self.dyes):]

        return G, color_map

    def graph_network(self, frame, include_dye=True):

        G, color_map = self.make_graph(frame, include_dye=include_dye)

        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_subplot(111)
        ax1.set_title(frame)
        nx.draw_networkx(G, with_labels=True, node_color=color_map)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        for pos in ['right','top','bottom','left']:
            plt.gca().spines[pos].set_visible(False)

        plt.savefig("network_" + str(frame) + ".png")
        plt.show()

    def network_valency(self, frame, include_dye=True):

        G, color_map = self.make_graph(frame, include_dye=include_dye)
        struct = list(nx.degree(G))

        poss = [struct[i][1] for i in range(len(G.nodes))]
         #       if struct[i] != 0]

        #print(list(nx.connected_components(G)))
        return sum(poss) / len(poss)

    def graph_clusters_all(self, frame, save_name=None):

        clusters = self.get_clusters_all(frame)

        sizes = [len(line) for line in clusters]

        #G, color_map = self.make_graph(frame, include_dye=include_dye)

        #clusters = list(nx.connected_components(G))

        #sizes = [len(list(cluster)) for cluster in clusters]

        occurences = [0 for _ in range(np.max(sizes) + 1)]

        for thing in sizes:
            occurences[thing] += 1
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title(frame)

        ax1.set_xlabel('size')
        ax1.set_ylabel('occurences')
        ax1.plot(list(range(0, np.max(sizes) + 1)), occurences)
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_clusters_from_data(self, all_sizes, save_name=None, combined=False):

        occurences = [0 for _ in range(np.max(all_sizes) + 2)]

        for thing in all_sizes:
            #print(thing, len(occurences))
            occurences[thing] += 1

        occurences = np.multiply(occurences, range(max(all_sizes) + 2))

        bins = np.zeros(100)

        for ind, value in enumerate(occurences):
            bins_index = int(np.floor((ind-1)/10))
            #print(ind, value)
            bins[bins_index] += value
        bins = np.divide(bins, np.sum(bins))
        occurences = np.divide(occurences, np.sum(occurences))

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")
        ax1.set_xlim([0,1])
        ax1.set_ylim([0, 1.2])
        ax1.set_xlabel('Fractional Cluster Size')
        ax1.set_ylabel('Weight Average Probability')
        #ax1.plot(list(range(1, np.max(all_sizes) + 2)), occurences[1:])
        ax1.plot(np.divide(list(range(10,1001,10)), 1000), bins)
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_clusters_frames(self, frames, save_name=None, include_dye=True):

        all_sizes = []

        for frame in frames:
            G, color_map = self.make_graph(frame, include_dye=include_dye)
            clusters = list(nx.connected_components(G))
            sizes = [len(list(cluster)) for cluster in clusters]
            #sizes = [len(line) for line in clusters]
            all_sizes.extend(sizes)
        self.graph_clusters_from_data(all_sizes, save_name=save_name)


    def get_all_aa_bonds(self, frame):

        g, cmap = self.make_graph(frame)
        edges = g.edges
        connections = []
        for edge in edges:
            first = edge[0]
            second = edge[1]
            connections.extend(self.get_specific_connection(self.all_chain_indices[first], self.all_chain_indices[second], frame))
        return connections

    def get_specific_connection(self, chain1, chain2, frame, cut=0.5):

        #print("checking new")
        cons = []
        real_frame = self.trajectory.read_frame(frame)
        positions1 = np.array([self.get_position(index1, real_frame) for mon1 in chain1 for index1 in mon1
                      if self.get_type(index1, real_frame) == "A"])
        positions2 = np.array([self.get_position(index2, real_frame) for mon2 in chain2 for index2 in mon2
                      if self.get_type(index2, real_frame) == "A"])
        indexes1 = [index1 for mon1 in chain1 for index1 in mon1 if self.get_type(index1, real_frame) == "A"]
        indexes2 = [index2 for mon2 in chain2 for index2 in mon2 if self.get_type(index2, real_frame) == "A"]

        if len(positions2) == 0 or len(positions1) == 0:
            return cons

        dist_array = distances.distance_array(positions1, positions2)
        s = dist_array.shape
        for i in range(s[0]):
            for j in range(s[1]):
                if dist_array[i][j] < cut:
                    cons.append([indexes1[i], indexes2[j]])
        return cons

    def get_broken_bonds(self, frame, cut=1):

        real_frame = self.trajectory.read_frame(frame)
        all_bonds = real_frame.bonds.group
        all_bond_types = real_frame.bonds.types

        broke = []

        for ind, bond in enumerate(all_bonds):
            if all_bond_types[real_frame.bonds.typeid[ind]][:8] == "go_morse":
                pos_1 = real_frame.particles.position[bond[0]]
                pos_2 = self.new_pos(real_frame.particles.position[bond[1]], pos_1, real_frame.configuration.box[:3], silence=True)
                dist = np.linalg.norm(np.subtract(pos_1, pos_2))
                if dist > cut:
                    broke.append([bond[0], bond[1]])
        return broke


    def graph_broken_bonds(self, frames, save_name="broken.png", strain_func=None):

        numbers = []
        for frame in frames:
            numbers.append(len(self.get_broken_bonds(frame)))
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")

        ax1.set_xlabel('time')
        ax1.set_ylabel('Broken Bonds')
        if strain_func is None:
            ax1.plot(frames, numbers)
        else:
            ax1.plot(strain_func(frames), numbers)
            ax1.set_xlabel("Strain ($L/L_{0}$)")
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)

        if strain_func is not None:
            self.write_file("broken.txt", strain_func(frames), numbers)
        plt.show()

    def convert_topology(self, frame, new_name):

        import hoomd
        hoomd.context.initialize("--mode=cpu")
        snap = hoomd.data.gsd_snapshot(self.gsd_name, frame=frame)

        bonds = self.get_all_aa_bonds(frame)

        k_orig = 150000
        k = k_orig / 2.479
        d = 140
        beta = np.sqrt(k / (2 * d))
        bond_name = "go_morse_0.38_" + str(d) + "_" + str(beta)

        print(snap.bonds)
        for bond in bonds:
            bond_number = snap.bonds.N + 1
            snap.bonds.resize(bond_number)
            snap.bonds.group[bond_number - 1] = bond
            snap.bonds.typeid[bond_number - 1] = 0

        import hoomd

        sys = hoomd.init.read_snapshot(snap)
        hoomd.dump.gsd(filename=new_name + ".gsd", period=None, group=hoomd.group.all(), overwrite=True)

    def get_protein_coordinates(self, frame, proteins=None):

        from MartiniProtein import MartiniProteinFromData

        real_frame = self.trajectory.read_frame(frame)
        if proteins is None:
            proteins=list(range(len(self.dyes)))
            for i in proteins:
                pos = np.array([self.get_position_min_image(index1, real_frame) for mon1 in self.dyes[i] for index1 in mon1])
                mass = np.array([real_frame.particles.mass[index1] for mon1 in self.dyes[i] for index1 in mon1])
                tipe = [self.get_type(index1, real_frame) for mon1 in self.dyes[i] for index1 in mon1]
                mp = MartiniProteinFromData(pos, tipe, mass)
                name = self.gsd_name[:-4] + "_protein" + str(i) + "_frame" + str(frame)
                mp.go_coordinate_file(name)

    def rmsd(self, protein_index, frame, reference_frame=0):


        if self.go_protein is None:
            real_reference_frame = self.trajectory.read_frame(reference_frame)
            ref_box = real_reference_frame.configuration.box[:3]
            reference = self.get_positions_chain([ind for ind in self.dyes[protein_index]], real_reference_frame, ref_box)
        else:
            reference = self.go_protein.position

        real_frame = self.trajectory.read_frame(frame)
        box = real_frame.configuration.box[:3]
        pos = self.get_positions_chain([ind for ind in self.dyes[protein_index]], real_frame, box)

        #if protein_index ==5 or protein_index == 8:
        #    print(pos-reference)
        return rms.rmsd(pos, reference, center=True, superposition=True)

    def broken_contacts(self, protein_index, frame):

        if self.go_protein is None:
            return False
        possible = self.go_protein.mp.contact_pairs
        real_frame = self.trajectory.read_frame(frame)
        box = real_frame.configuration.box[:3]
        pos = self.get_positions_chain([ind for ind in self.dyes[protein_index]], real_frame, box)
        broken = []

        for contact in possible:
            ind1 = self.go_protein.mp.backbone_indices.index(contact[0])
            ind2 = self.go_protein.mp.backbone_indices.index(contact[1])
            c6 = contact[3]
            c12 = contact[4]

            sigma = (c6 / c12) ** (-1 / 6)
            dist = np.linalg.norm(np.subtract(pos[ind1], pos[ind2]))
            if dist > 1.1:
                broken.append(1)
                #print(ind1, ind2, 1)
                #print(contact, sigma * 2**(1/6), dist)
            else:
                broken.append(0)
                #print(ind1, ind2, 0)
                #print(contact, sigma * 2 **(1/6), dist)
        return broken

    def average_broken_contacts(self, frame):

        broken = []

        for i in range(len(self.dyes)):
            broken.extend(self.broken_contacts(i, frame))

        #print(frame, np.sum(broken) / len(broken))
        return np.sum(broken) / len(broken)


    def dump_protein_xyz(self, protein_index, frame, name=None, centered=True):

        real_frame = self.trajectory.read_frame(frame)
        box = real_frame.configuration.box[:3]
        pos = self.get_positions_chain([ind for ind in self.dyes[protein_index]], real_frame, box)

        center = np.average(pos, axis=0)
        pos = np.subtract(pos, center)
        if name is None:
            name = "protein_" + str(protein_index) + ".xyz"
        f = open(name, "w")
        f.write(str(len(pos)) + "\n\n")
        for i in range(len(pos)):
            if i ==0:
                tipe = "N"
            elif i == len(pos) - 1:
                tipe = "C"
            else:
                tipe = "BB"
            f.write(tipe + " " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n")



    def graph_rmsd_dist(self, frame, save_name=None):

        numbers = []
        for i in range(len(self.dyes)):
            numbers.append(self.rmsd(i, frame, reference_frame=0))

        bin_size=.2
        hist, bin_edges = np.histogram(numbers, bins=int(np.max(numbers)/bin_size))


        bins = [(bin_edges[i] + bin_edges[i+1]) /2 for i in range(len(bin_edges) -1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")

        ax1.set_xlabel('RMSD')
        ax1.set_ylabel('Probability')


        ax1.plot(bins, hist)
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
        else:
            save_name = str(frame) + "_rmsd_dist.png"

        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_average_rmsd(self, frames, save_name="rmsd_ave.png", strain_func=None):

        numbers = []
        for frame in frames:
            numbers.append(np.average([self.rmsd(i, frame, reference_frame=0) for i in range(len(self.dyes))]))
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")

        ax1.set_xlabel('time')
        ax1.set_ylabel('Average RMSD')
        if strain_func is None:
            ax1.plot(frames, numbers)
        else:
            ax1.plot(strain_func(frames), numbers)
            ax1.set_xlabel("Strain ($L/L_{0}$)")
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        if strain_func is not None:
            self.write_file("rmsd.txt", strain_func(frames), numbers)
        plt.show()

    def write_file(self, name, x, y):

        f = open(name, "w")

        for ind, x_thing in enumerate(x):
            f.write(str(x_thing) + " " + str(y[ind]) + "\n")


    def get_positions_chain(self, indices, frame, box):

        positions = []
        for number, index in enumerate(indices):
            position = frame.particles.position[index]
            if number != 0:
                position = self.new_pos(position, frame.particles.position[indices[number-1]], box)
            positions.append(position)

        return positions


    def new_pos(self, pos, ref, box,  silence=False):

        do = False
        if la.norm(np.subtract(pos, ref)) > 10:
        #    print("start")
        #    print(pos, ref)
            do = True
        x_dist = pos[0] - ref[0]
        if x_dist > box[0] / 2:
            pos[0] = pos[0] - box[0]
        elif x_dist < - box[0]/2:
            pos[0] = pos[0] + box[0]

        y_dist = pos[1] - ref[1]
        if y_dist > box[1] / 2:
            pos[1] = pos[1] - box[1]
        elif y_dist < - box[1]/2:
            pos[1] = pos[1] + box[1]

        z_dist = pos[2] - ref[2]
        if z_dist > box[2] / 2:
            pos[2] = pos[2] - box[2]
        elif z_dist < - box[2]/2:
            pos[2] = pos[2] + box[2]
        if do:
         #   print(pos, ref)
         #   print(la.norm(np.subtract(pos, ref)))
            if la.norm(np.subtract(pos,ref))>10 and not silence:
                print("aaaaaaaa")
         #   print("done")
        return pos


