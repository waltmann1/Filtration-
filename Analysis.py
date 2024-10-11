from __future__ import division
import gsd.hoomd
import gsd.fl
import numpy as np
import numpy.linalg as la
import os.path
import networkx as nx
import copy as cp
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})
import hoomd.data as hbox
from mpl_toolkits.mplot3d import Axes3D
import math as m
import MDAnalysis as mda
from MDAnalysis.analysis import distances


class Analysis(object):

    def __init__(self, gsd_name, map_name):
        f = gsd.fl.open(name=gsd_name, mode='rb', application='', schema='hoomd',
                        schema_version=[1, 0])
        self.trajectory = gsd.hoomd.HOOMDTrajectory(f)
        self.dyes = None
        self.rando_chain_sequences = None
        self.rando_chain_indices = None
        self.charged_chain_indices = None
        self.charged_chain_sequences = None
        self.qpi = None
        self.qmi = None
        self.read_map(map_name)
        self.frames = []
        self.box = self.trajectory.read_frame(0).configuration.box[:3]
        self.gsd_name = gsd_name

    def get_box(self, frame):
        self.box = self.trajectory.read_frame(frame).configuration.box[:3]

    def read_map(self, map_name):

        f = open(map_name, 'r')

        dyes = []
        rando_chain_sequences = []
        rando_chain_indices = []
        charged_chain_indices = []
        charged_chain_sequences = []

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
                rando_chain_indices.append(chain)
                rando_chain_sequences.append(sequence)
            elif s[0] == "charged_polymer":
                i = 2
                sequence = [s[1]]
                mon = []
                chain = []
                while i < len(s):
                    if not s[i].isalpha():
                        mon.append(int(float(s[i])))
                    else:
                        chain.append(mon)
                        mon = []
                        sequence.append(s[i])
                    i += 1
                chain.append(mon)
                charged_chain_indices.append(chain)
                charged_chain_sequences.append(sequence)
            elif s[0] == "qpi":
                qpi = [int(float(s[i])) for i in range(1, len(s))]
            elif s[0] == "qmi":
                qmi = [int(float(s[i])) for i in range(1, len(s))]

        self.dyes = dyes
        self.rando_chain_sequences = rando_chain_sequences
        self.rando_chain_indices = rando_chain_indices
        self.charged_chain_indices = charged_chain_indices
        self.charged_chain_sequences = charged_chain_sequences
        self.all_chain_indices = rando_chain_indices + charged_chain_indices
        self.all_chain_sequences = rando_chain_sequences + charged_chain_sequences
        self.qpi = qpi
        self.qmi = qmi

        self.charged_mers = ['SPMMA', 'QVP', 'MAETMA', 'SPMMA2']

        self.hydrophobic_mers = ['EHMA', 'VP']

        self.hydrophilic_mers = ['PEGMEMA']

    def chain_rgs_backbone(self, frame, positives=False, negatives=True, max_hydro=2.0, min_hydro=-1.0, in_cluster=None):

        self.get_box(frame)
        frame_number =frame
        frame = self.trajectory.read_frame(frame)
        rgs = []
        chain_list = self.all_chain_indices
        indices = []
        if positives and not negatives:
            chain_list = chain_list[len(self.rando_chain_indices):]
            indices = list(range(len(self.rando_chain_indices), len(self.all_chain_indices)))
        elif negatives and not positives:
            chain_list = chain_list[:len(self.rando_chain_indices)]
            indices = list(range(len(self.rando_chain_indices)))
        elif negatives and positives:
            chain_list = chain_list
            indices = list(range(len(self.all_chain_indices)))
        if in_cluster is not None:
            list_in_cluster = [self.in_cluster(i, frame_number) for i in indices]
            for truth_index in range(len(list_in_cluster)-1, -1, -1):
                if list_in_cluster[truth_index] != in_cluster:
                    indices.remove(indices[truth_index])
                    chain_list.remove(chain_list[truth_index])
        for i, chain in enumerate(chain_list):
            index = indices[i]
            percent_hydro = self.get_percent_hydrophobic(index)
            if percent_hydro <= max_hydro and percent_hydro > min_hydro:
                backbone = [mon[0] for mon in chain]
                positions = self.get_positions_chain(backbone, frame)
                com = np.average(positions, axis=0)
                squares = [self.distance(com, pos)**2 for pos in positions]
                rgs.append(np.sqrt(np.sum(squares)/len(backbone)))
        return rgs

    def get_positions_chain(self, indices, frame):

        positions = []
        for number, index in enumerate(indices):
            position = frame.particles.position[index]
            if number != 0:
                position = self.new_pos(position, frame.particles.position[indices[number-1]])
            positions.append(position)

        return positions

    def rdf_by_chain(self, frame, cut=1000):

        self.get_box(frame)
        frame_number = frame
        frame = self.trajectory.read_frame(frame)
        types = [self.get_type(index, frame) for index in range(len(frame.particles.position))]
        types = np.unique(types)
        types = list(types)
        rdf_distances = [[] for _ in range(len(types))]
        for chain in self.all_chain_indices:
            #print(chain, chain[0])
            indices = [i for mon in chain for i in mon]
            positions = self.get_positions_chain(indices, frame)
            com = np.average(positions, axis=0)
            distances = la.norm(np.subtract(positions, com), axis=1)

            for ind, distance in enumerate(distances):
                self.total += 1
                if np.any(distances > cut):
                    self.counted += 1
                    rdf_distances[types.index(self.get_type(indices[ind], frame))].append(distance)

        return types, rdf_distances

    def graph_micelle_rdf(self, frames, cut=1000):

        self.total = 0
        self.counted = 0
        types, rdf_distances = self.rdf_by_chain(frames[0], cut=cut)
        #print(len(types), len(rdf_distances))
        #print(len(types[0]), len(rdf_distances[0]))
        #print(types)
        #quit()
        for frame in frames[1:]:
            new_types, new_r_d = self.rdf_by_chain(frame, cut=cut)
            for i in range(len(new_types)):
                rdf_distances[types.index(new_types[i])].extend(rdf_distances[i])

        bins = [1 * i for i in range(22)]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('R (nm)')
        ax1.set_ylabel('P(R)')

        colors = ["#aaaa00", "#0000ff", "#00ff7f"]

        names = ["NonPolar", "Polar", "Negatively Charged"]


        keepers = rdf_distances[:2]
        keepers.append(rdf_distances[3])

        for ind, rdf in enumerate(keepers):
            hist, bin_edges = np.histogram(rdf, bins=bins)
            #if names[ind] == "Polar":
            #    hist = np.multiply(hist, 1.5)
            hist = np.divide(hist, np.sum(hist))
            #print(names[ind], hist)
            ax1.plot(np.multiply(bin_edges[:-1], .5), hist, label=names[ind], color=colors[ind])

        print(self.total, self.counted, self.counted/self.total)
        #plt.legend(fontsize=14, loc="upper right", ncol=1, bbox_to_anchor=(1.1, 1.0))
        plt.savefig("micelle_hist.png", bbox_inches='tight', pad_inches=.2)
        plt.show()


    def new_pos(self, pos, ref, silence=False):

        do = False
        if la.norm(np.subtract(pos, ref)) > 10:
        #    print("start")
        #    print(pos, ref)
            do = True
        x_dist = pos[0] - ref[0]
        if x_dist > self.box[0] / 2:
            pos[0] = pos[0] - self.box[0]
        elif x_dist < - self.box[0]/2:
            pos[0] = pos[0] + self.box[0]

        y_dist = pos[1] - ref[1]
        if y_dist > self.box[1] / 2:
            pos[1] = pos[1] - self.box[1]
        elif y_dist < - self.box[1]/2:
            pos[1] = pos[1] + self.box[1]

        z_dist = pos[2] - ref[2]
        if z_dist > self.box[2] / 2:
            pos[2] = pos[2] - self.box[2]
        elif z_dist < - self.box[2]/2:
            pos[2] = pos[2] + self.box[2]
        if do:
         #   print(pos, ref)
         #   print(la.norm(np.subtract(pos, ref)))
            if la.norm(np.subtract(pos,ref))>10 and not silence:
                print("aaaaaaaa")
         #   print("done")
        return pos

    def distance(self, one, two):
        #if la.norm(np.subtract(one, two)) > 10:
         #   print(one, two)
        #return np.abs(np.subtract(one[2], two[]))
        #new_pos = self.new_pos(one, two, silence=True)
        return la.norm(np.subtract(one, two))


    def graph_rg_distribution(self, frame):

        dists = self.chain_rgs_backbone(frame)
        rdf_hist, rbe = np.histogram(dists, bins=10)
        bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Rg')
        average = np.average(dists)
        ax1.set_ylabel('')
        ax1.plot(bin_middles, rdf_hist, label='<Rg = ' + str(average)+ '>')
        plt.legend()
        plt.show()

    def graph_combined_rg_distribution(self, frames, save_name=None, positives=False, negatives=True,
                                       max_hydro=2.0, min_hydro=-1.0, in_cluster=None):

        dists = []
        for frame in frames:
            dists += self.chain_rgs_backbone(frame, positives=positives, negatives=negatives,
                                             max_hydro=max_hydro, min_hydro=min_hydro, in_cluster=in_cluster)
        dists = [dist for dist in dists if dist < 18]
        rdf_hist, rbe = np.histogram(dists, bins=10)
        bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Rg')
        average = np.average(dists)
        ax1.set_ylabel('')
        ax1.plot(bin_middles, rdf_hist, label='<Rg = ' + str(average) + '>')
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_clusters_all(self, frame, save_name=None):

        clusters = self.get_clusters_all(frame)

        sizes = [len(line) for line in clusters]

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


    def graph_clusters_frames(self, frames, save_name=None):

        all_sizes = []

        for frame in frames:
            clusters = self.get_clusters_all(frame)
            sizes = [len(line) for line in clusters]
            all_sizes.extend(sizes)
        self.graph_clusters_from_data(all_sizes, save_name=save_name)

    def graph_clusters_from_data(self, all_sizes, save_name=None, combined=False):

        occurences = [0 for _ in range(np.max(all_sizes) + 2)]

        for thing in all_sizes:
            #print(thing, len(occurences))
            occurences[thing] += 1

        occurences = np.multiply(occurences, range(max(all_sizes) + 2))

        bins = np.zeros(12)

        for ind, value in enumerate(occurences):
            bins_index = int(np.floor((ind-1)/5))
            print(ind, value)
            bins[bins_index] += value
        bins = np.divide(bins, np.sum(bins))
        occurences = np.divide(occurences, np.sum(occurences))

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("")
        ax1.set_xlim([0,1])
        ax1.set_ylim([0, .5])
        ax1.set_xlabel('Fractional Cluster Size')
        ax1.set_ylabel('Weight Average Probability')
        #ax1.plot(list(range(1, np.max(all_sizes) + 2)), occurences[1:])
        ax1.plot(np.divide(list(range(5,61,5)), 60), bins)
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()


    def graph_clusters(self, frame, save_name=None):

        clusters = self.get_clusters(frame)

        sizes = [len(line) for line in clusters]

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

    def graph_network(self, frame, include_pos=True):

        G, color_map = self.make_graph(frame, include_pos=include_pos)

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

    def network_valency_positives(self, frame):

        G, color_map = self.make_graph(frame, include_pos=True)
        struct = nx.degree(G)
        poss = [struct[i] for i in range(len(self.rando_chain_indices), len(self.all_chain_indices))
                if struct[i] != 0]

        return sum(poss) / len(poss)

    def network_valency_negatives(self, frame):

        G, color_map = self.make_graph(frame, include_pos=True)
        struct = nx.degree(G)
        #negs = [struct[i] for i in range(0, len(self.rando_chain_indices))
        #        if struct[i] != 0]
        negs = [struct[i] for i in range(0, len(self.rando_chain_indices))]
        return sum(negs) / len(negs)

    def graph_network_valencies_frames(self, frames, save_name=None):

        positives = []
        negatives = []
        for frame in frames:
            positives.append(self.network_valency_positives(frame))
            negatives.append(self.network_valency_negatives(frame))

        fig = plt.figure()

        ax1 = fig.add_subplot(111)


        ax1.set_xlabel('Frame')
        ax1.set_ylabel('Average Valency of Connected Chains')
        ax1.plot(frames, negatives, label='Negatives', color='blue')
        ax1.plot(frames, positives, label='Positives', color='r')
        print(np.average(positives))
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_number_of_connections(self, frames, save_name=None):

        connect = []
        neg_neg_list=[]
        pos_pos_list=[]
        neg_pos_list=[]
        neg_number = len(self.rando_chain_sequences)
        for frame in frames:
            neg_neg=0
            pos_pos=0
            neg_pos=0
            G, color_map = self.make_graph(frame, include_pos=True)
            connect.append(G.number_of_edges())
            for edge in G.edges:
                if edge[0] < neg_number and edge[1] < neg_number:
                    neg_neg += 1
                elif (edge[0] > neg_number and edge[1] > neg_number):
                    pos_pos += 1
                else:
                    neg_pos += 1
            neg_neg_list.append(neg_neg)
            neg_pos_list.append(neg_pos)
            pos_pos_list.append(pos_pos)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)


        ax1.set_xlabel('Frame')
        ax1.set_ylabel('Connections')
        ax1.set_ylim([0, 100])
        ax1.plot(frames, connect, label='Total', color='k')
        ax1.plot(frames, neg_neg_list, label='Negative-Negative', color='blue')
        ax1.plot(frames, neg_pos_list, label='Positive-Negative', color='r')
        ax1.plot(frames, pos_pos_list, label='Positive-Positive', color='pink')
        plt.legend()
        print(np.average(connect))
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def make_graph(self, frame, include_pos=True):

        G = nx.Graph()
        color_map = []
        total_number = len(self.all_chain_indices)
        neg_number = len(self.rando_chain_indices)
        pos_number = total_number - neg_number
        for i in range(total_number):
            G.add_node(i)
            if i < neg_number:
                color_map.append('blue')
            else:
                color_map.append('r')

        clusters = self.get_all_connections_from_clusters(frame)

        for i in range(neg_number + include_pos * pos_number):
            for thing in clusters[i]:
                if thing < neg_number + include_pos * pos_number and i != thing:
                    G.add_edge(i, thing)
        return G, color_map

    def get_percent_charged(self, index):

        count = 0
        sequence = self.all_chain_sequences[index]
        for thing in sequence:
            if thing in self.charged_mers:
                count += 1
        return count / len(sequence)

    def get_percent_hydrophobic(self, index):

        count = 0
        sequence = self.all_chain_sequences[index]
        for thing in sequence:
            if thing in self.hydrophobic_mers:
                count += 1
        return count / len(sequence)

    def get_percent_hydrophilic(self, index):

        count = 0
        sequence = self.all_chain_sequences[index]
        for thing in sequence:
            if thing in self.hydrophilic_mers:
                count += 1
        return count / len(sequence)

    def in_cluster(self, index, frame):

        clusters = self.get_clusters_all(frame)
        #clusters = self.get_clusters(frame)

        for cluster in clusters:
            if index in cluster:
                return len(cluster) == 1
        print("Warning: index does not exist")

    def percent_charged_in_clusters(self, frame, negatives=True, positives=False):

        clusters = self.get_clusters_all(frame)

        in_cluster = []
        not_in_cluster = []
        for cluster in clusters:
            if len(cluster) > 1:
                for index in cluster:
                    if index >= len(self.rando_chain_sequences) and positives:
                        in_cluster.append(self.get_percent_charged(index))
                    if index < len(self.rando_chain_sequences) and negatives:
                        in_cluster.append(self.get_percent_charged(index))
            else:
                index = cluster[0]
                if index >= len(self.rando_chain_sequences) and positives:
                    not_in_cluster.append(self.get_percent_charged(index))
                if index < len(self.rando_chain_sequences) and negatives:
                    not_in_cluster.append(self.get_percent_charged(index))

        i_hist, ibe = np.histogram(in_cluster, bins=10)
        bin_middles = [(ibe[i] + ibe[i + 1]) / 2 for i in range(len(ibe) - 1)]

        ni_hist, nibe = np.histogram(not_in_cluster, bins=10)
        nbin_middles = [(nibe[i] + nibe[i + 1]) / 2 for i in range(len(nibe) - 1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Fq')
        ax1.set_ylabel('')
        ax1.plot(bin_middles, i_hist, label='in clusters')
        ax1.plot(nbin_middles, ni_hist, label="not in clusters")
        plt.legend()
        plt.show()

    def percent_hydrophobic_in_clusters(self, frame, negatives=True, positives=False):

        clusters = self.get_clusters_all(frame)

        in_cluster = []
        not_in_cluster = []
        for cluster in clusters:
            if len(cluster) > 1:
                for index in cluster:
                    if index >= len(self.rando_chain_sequences) and positives:
                        in_cluster.append(self.get_percent_hydrophobic(index))
                    if index < len(self.rando_chain_sequences) and negatives:
                        in_cluster.append(self.get_percent_hydrophobic(index))
            else:
                index = cluster[0]
                if index >= len(self.rando_chain_sequences) and positives:
                    not_in_cluster.append(self.get_percent_hydrophobic(index))
                if index < len(self.rando_chain_sequences) and negatives:
                    not_in_cluster.append(self.get_percent_hydrophobic(index))

        i_hist, ibe = np.histogram(in_cluster, bins=10)
        bin_middles = [(ibe[i] + ibe[i + 1]) / 2 for i in range(len(ibe) - 1)]

        ni_hist, nibe = np.histogram(not_in_cluster, bins=10)
        nbin_middles = [(nibe[i] + nibe[i + 1]) / 2 for i in range(len(nibe) - 1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Fq')
        ax1.set_ylabel('')
        ax1.plot(bin_middles, i_hist, label='in clusters')
        ax1.plot(nbin_middles, ni_hist, label="not in clusters")
        plt.legend()
        plt.show()

    def percent_hydrophilic_in_clusters(self, frame, negatives=True, positives=False):

        clusters = self.get_clusters_all(frame)

        in_cluster = []
        not_in_cluster = []
        for cluster in clusters:
            if len(cluster) > 1:
                for index in cluster:
                    if index >= len(self.rando_chain_sequences) and positives:
                        in_cluster.append(self.get_percent_hydrophilic(index))
                    if index < len(self.rando_chain_sequences) and negatives:
                        in_cluster.append(self.get_percent_hydrophilic(index))
            else:
                index = cluster[0]
                if index >= len(self.rando_chain_sequences) and positives:
                    not_in_cluster.append(self.get_percent_hydrophilic(index))
                if index < len(self.rando_chain_sequences) and negatives:
                    not_in_cluster.append(self.get_percent_hydrophilic(index))

        i_hist, ibe = np.histogram(in_cluster, bins=10)
        bin_middles = [(ibe[i] + ibe[i + 1]) / 2 for i in range(len(ibe) - 1)]

        ni_hist, nibe = np.histogram(not_in_cluster, bins=10)
        nbin_middles = [(nibe[i] + nibe[i + 1]) / 2 for i in range(len(nibe) - 1)]

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Fq')
        ax1.set_ylabel('')
        ax1.plot(bin_middles, i_hist, label='in clusters')
        ax1.plot(nbin_middles, ni_hist, label="not in clusters")
        plt.legend()
        plt.show()

    def get_clusters(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters.txt")

        not_in_cluster = list(range(0, len(self.rando_chain_indices)))

        clusters = []

        current_cluster =[]

        while len(not_in_cluster) != 0:
            print("not in cluster", not_in_cluster)
            self.look_for_connections(not_in_cluster, current_cluster, not_in_cluster[0], frame)
            clusters.append(current_cluster)
            current_cluster = []

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters.txt", "w")
        for cluster in clusters:
            string = ""
            for thing in cluster:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return clusters

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

    def get_clusters_charged(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters3.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters3.txt")

        #not_in_cluster = list(range(0, len(self.all_chain_indices)))

        clusters = []

        current_cluster = []

        for index1, chain_1 in enumerate(self.all_chain_indices):
            current_cluster.append(index1)
            for index2, chain_2 in enumerate(self.all_chain_indices):
                if self.look_for_connection_charged(chain_1, chain_2, frame):
                    current_cluster.append(index2)
                    print("hit " + str(index1) + " " + str(index2))
            clusters.append(current_cluster)
            current_cluster = []
            #print("clusters ", clusters)

        #while len(not_in_cluster) != 0:
        #    print("not in cluster", not_in_cluster)
        #    self.look_for_connections_charged(self.all_chain_indices, frame)
        #    clusters.append(current_cluster)
         #   current_cluster = []

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters3.txt", "w")
        for cluster in clusters:
            string = ""
            for thing in cluster:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return clusters

    def get_all_connections_from_clusters(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters4.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_clusters4.txt")

        #not_in_cluster = list(range(0, len(self.all_chain_indices)))

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

    def look_for_connections(self, not_in_cluster,current_cluster,current_index, frame):
        print(current_index)
        current_cluster.append(current_index)
        not_in_cluster.remove(current_index)
        for chain_index in not_in_cluster:
            if self.look_for_connection(self.rando_chain_indices[current_index], self.rando_chain_indices[chain_index], frame):
                print("hit", current_index, chain_index)
                self.look_for_connections(not_in_cluster, current_cluster, chain_index, frame)

    def look_for_connections_all(self, not_in_cluster,current_cluster,current_index, frame):
        print(current_index)
        current_cluster.append(current_index)
        not_in_cluster.remove(current_index)
        for chain_index in not_in_cluster:
            if self.look_for_connection(self.all_chain_indices[current_index], self.all_chain_indices[chain_index], frame):
                print("hit", current_index, chain_index)
                self.look_for_connections_all(not_in_cluster, current_cluster, chain_index, frame)

    def look_for_connections_charged(self, chains, frame):
        for chain1 in chains:
            for chain2 in chains:
                return self.look_for_connection_charged(chain1, chain2, frame)

    def look_for_connection(self, chain1, chain2, frame, cut=1.122):

        print("checking new")
        #for mon1 in chain1:
        #    for index1 in mon1:
        #        if self.get_type(index1, frame)[0] == "B":
        #            for mon2 in chain2:
        #                for index2 in mon2:
        #                    if self.get_type(index2, frame)[0] == "B":
        #                        if self.distance(self.get_position(index1, frame), self.get_position(index2, frame)) \
        #                                < cut:
        #                            return True
        positions1 = np.array([self.get_position(index1, frame) for mon1 in chain1 for index1 in mon1
                      if self.get_type(index1, frame) == "B"])
        positions2 = np.array([self.get_position(index2, frame) for mon2 in chain2 for index2 in mon2
                      if self.get_type(index2, frame) == "B"])
        dist_array = distances.distance_array(positions1, positions2)
        if np.any(dist_array < cut):
            return True



        return False

    def look_for_connection_charged(self, chain1, chain2, frame, cut=5):

        #print("checking charged")
        lb = .7
        for mon1 in chain1:
            for index1 in mon1:
                energy = 0
                tipe1 = self.get_type(index1, frame)
                if self.get_type(index1, frame)[0] == "Q" and self.get_type(index1, frame)[-1] != 'i':
                    if tipe1 == "QP":
                        return False
                    for mon2 in chain2:
                        for index2 in mon2:
                            if self.get_type(index2, frame)[0] == "Q" and self.get_type(index2, frame)[-1] != 'i':
                                if self.get_type(index2, frame) == "QM":
                                    return False
                                if self.distance(self.get_position(index1, frame), self.get_position(index2, frame)) < 5* lb:
                                    energy += lb / self.distance(self.get_position(index1, frame), self.get_position(index2, frame))
                    #print(energy)
                    if energy > 2:
                        print("index " + str(index1))
                        return True
        return False

    def get_monomer_connections(self, chain1, chain2, frame, cut=1.122):

        chain1_monomers = []
        chain2_monomers = []

        chain1_indices = np.array([mon1 for mon1 in chain1 for index1 in mon1
                               if self.get_type(index1, frame) == "B"])

        chain2_indices = np.array([mon2 for mon2 in chain2 for index2 in mon2
                               if self.get_type(index2, frame) == "B"])

        positions1 = np.array([self.get_position(index1, frame) for mon1 in chain1 for index1 in mon1
                               if self.get_type(index1, frame) == "B"])
        positions2 = np.array([self.get_position(index2, frame) for mon2 in chain2 for index2 in mon2
                               if self.get_type(index2, frame) == "B"])
        dist_array = distances.distance_array(positions1, positions2)


        for ind_row, dist_row in enumerate(dist_array):
            if np.any(dist_row < cut) and chain1_indices[ind_row] not in chain1_monomers:
                chain1_monomers.append(chain1_indices[ind_row])
            for ind, thing in enumerate(dist_row):
                if thing < cut and chain2_indices[ind] not in chain2_monomers:
                    chain2_monomers.append(chain2_indices[ind])

        print(chain1_monomers, chain2_monomers)
        return chain1_monomers, chain2_monomers

    def get_clustered_monomer_sequences(self, frame):

        clusters = self.get_clusters_all(frame)
        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences.txt", seq=True)

        sequences = []
        for cluster in clusters:
            print(cluster)
            if len(cluster) > 1:
                for ind, chain in enumerate(cluster):
                    for chain2 in cluster[ind+1:]:
                        use_chain = self.all_chain_indices[chain]
                        chain2 = self.all_chain_indices[chain2]
                        indexes1, indexes2 = self.get_monomer_connections(use_chain, chain2, frame)
                        print(len(indexes1), len(indexes2))
                        for index in indexes1:
                            mon_index = self.get_monomer_index_from_tag(index, use_chain)
                            sequences.append(self.get_sequence_around_monomer(mon_index,
                                                                              self.all_chain_indices.index(use_chain)))
                        for index in indexes2:
                            mon_index = self.get_monomer_index_from_tag(index, chain2)
                            sequences.append(self.get_sequence_around_monomer(mon_index,
                                                                              self.all_chain_indices.index(chain2)))
        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences.txt", "w")
        for sequence in sequences:
            string = ""
            for thing in sequence:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return sequences

    def get_pos_neg_sequences(self, frame):

        clusters = self.get_clusters_all(frame)
        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences2.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences2.txt", seq=True)

        sequences = []
        for cluster in clusters:
            print(cluster)
            neg_num = len(self.rando_chain_indices)
            if len(cluster) > 1:
                for ind, chain in enumerate(cluster):
                    for chain2 in cluster[ind+1:]:
                        if (chain2 >= neg_num and chain < neg_num) or (chain2 < neg_num and chain >= neg_num):
                            use_chain = self.all_chain_indices[chain]
                            chain2 = self.all_chain_indices[chain2]
                            indexes1, indexes2 = self.get_monomer_connections(use_chain, chain2, frame)
                            print(len(indexes1), len(indexes2))
                            for index in indexes1:
                                mon_index = self.get_monomer_index_from_tag(index, use_chain)
                                sequences.append(self.get_sequence_around_monomer(mon_index,
                                                                              self.all_chain_indices.index(use_chain)))
                            for index in indexes2:
                                mon_index = self.get_monomer_index_from_tag(index, chain2)
                                sequences.append(self.get_sequence_around_monomer(mon_index,
                                                                              self.all_chain_indices.index(chain2)))
        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_sequences2.txt", "w")
        for sequence in sequences:
            string = ""
            for thing in sequence:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        #f.close()
        return sequences

    def percent_hydrophobic_in_sequences(self, frame, save_name=None):

        #sequences = self.get_clustered_monomer_sequences(frame)

        sequences = []
        print("type of frame", type(frame))
        if True:
            for the_frame in frame:
                print(frame)
                x = self.get_pos_neg_sequences(the_frame)
                for thing in x:
                    sequences.append(thing)
        else:
            sequences = self.get_pos_neg_sequences(frame)

        num_seq = 0
        hydro = 0

        spmma = [0,0,0,0,0]
        ehma = [0,0,0,0,0]
        pegmema = [0,0,0,0,0]
        empty = [0, 0, 0, 0, 0]
        for sequence in sequences:
            if not ('BC' in sequence or 'QVP' in sequence or 'MAETMA' in sequence or 'MMA' in sequence) and len(sequence) > 0:
                num_seq += 1
                #if sequence.count("SPMMA") == 0:
                 #   hydro += 1
                for ind, mon in enumerate(sequence):
                    if mon[:5] == "SPMMA":
                        spmma[ind] += 1
                    if mon == "EHMA":
                        ehma[ind] += 1
                    if mon == "PEGMEMA":
                        pegmema[ind] += 1
                    if mon == "empty":
                        empty[ind] += 1

        spmma = np.divide(spmma, num_seq)
        ehma = np.divide(ehma, num_seq)
        pegmema = np.divide(pegmema, num_seq)
        #empty = np.divide(empty, num_seq)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Position')
        ax1.set_ylabel('Probability')
        ax1.plot(range(-2,3), spmma, label='SPMMA', color='blue')
        ax1.plot(range(-2,3), [.08] * 5, ls='dashed', color='blue')
        ax1.plot(range(-2,3), ehma, label='EHMA', color='purple')
        ax1.plot(range(-2,3), [.46] * 5, ls='dashed', color='purple')
        ax1.plot(range(-2,3), pegmema, label='PEGMEMA', color='green')
        ax1.plot(range(-2,3), [.46] * 5, ls='dashed', color='green')
        #ax1.plot(range(5), empty, label='Empty - chain terminated', color='k')
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()



    def percent_hydrophobic_in_dye_sequences(self, frames, save_name=None, negatives=True, positive=False):

        sequences = []
        for the_frame in frames:
            x = self.condensed_dye_sequences(the_frame)
            for thing in x:
                sequences.append(thing)
        num_f=len(frames)

        num_seq = 0
        hydro = 0

        spmma = [0,0,0,0,0]
        ehma = [0,0,0,0,0]
        pegmema = [0,0,0,0,0]
        empty = [0, 0, 0, 0, 0]
        pos = [0, 0, 0, 0, 0]
        other = [0, 0, 0, 0, 0]
        for sequence in sequences:
            if not ('VP' in sequence or 'QVP' in sequence or 'MAETMA' in sequence or 'MMA' in sequence) and len(sequence) > 0:
                num_seq += 1
                #if sequence.count("SPMMA") == 0:
                 #   hydro += 1
                for ind, mon in enumerate(sequence):
                    if mon[:5] == "SPMMA":
                        spmma[ind] += 1
                    if mon == "EHMA":
                        ehma[ind] += 1
                    if mon == "PEGMEMA":
                        pegmema[ind] += 1
                    if mon == "empty":
                        empty[ind] += 1

        print(num_seq, "in_negative", num_seq/num_f)
        spmma = np.divide(spmma, num_seq)
        ehma = np.divide(ehma, num_seq)
        pegmema = np.divide(pegmema, num_seq)
        #empty = np.divide(empty, num_seq)
        num_seq = 0
        for sequence in sequences:
            if ('VP' in sequence or 'QVP' in sequence or 'MAETMA' in sequence or 'MMA' in sequence) and len(sequence) > 0:
                num_seq += 1
                for ind, mon in enumerate(sequence):
                    if mon == "QVP" or mon == "MAETMA":
                        pos[ind] += 1
                    else:
                        other[ind] += 1
        pos = np.divide(pos, num_seq)
        other = np.divide(other, num_seq)
        print(num_seq, "in_positiive", num_seq/num_f)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Position')
        ax1.set_ylabel('Probability')
        if negatives:
            ax1.plot(range(-2,3), spmma, label='SPMMA', color='blue')
            ax1.plot(range(-2,3), [self.fraction_of_random("SPMMA2")] * 5, ls='dashed', color='blue')
            ax1.plot(range(-2,3), ehma, label='EHMA', color='purple')
            ax1.plot(range(-2,3), [self.fraction_of_random("EHMA")] * 5, ls='dashed', color='purple')
            ax1.plot(range(-2,3), pegmema, label='PEGMEMA', color='green')
            ax1.plot(range(-2,3), [self.fraction_of_random("PEGMEMA")] * 5, ls='dashed', color='green')
        if positive:
            ax1.plot(range(-2,3), pos, label=positive, color='red')
            ax1.plot(range(-2,3), [self.fraction_of_positive(positive)] * 5, ls='dashed', color='red')
            ax1.plot(range(-2,3), other, label='OTHER', color='yellow')
            ax1.plot(range(-2,3), [1 - self.fraction_of_positive(positive)] * 5, ls='dashed', color='yellow')

        #ax1.plot(range(5), empty, label='Empty - chain terminated', color='k')
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()


    def get_monomer_index_from_tag(self, tag, chain):

        return chain.index(tag)

    def fraction_of_random(self, monomer):
        total = np.sum([len(i) for i in self.rando_chain_indices])
        count = 0
        for chain in self.rando_chain_sequences:
            for mon in chain:
                if mon == monomer:
                    count += 1
        return count/total

    def fraction_of_positive(self, monomer):
        total = np.sum([len(i) for i in self.charged_chain_indices])
        count = 0
        for chain in self.charged_chain_sequences:
            for mon in chain:
                if mon == monomer:
                    count += 1
        return count/total

    def get_sequence_around_monomer(self, monomer_index, chain_index):

        number_around =2
        pad = ["empty"] * number_around
        chain = self.all_chain_sequences[chain_index]
        if monomer_index < number_around:
            chain = pad + chain
            #print(chain[:5])
            sequence = chain[monomer_index: monomer_index + 2*number_around + 1]
        else:
            if (len(chain) - 1) - monomer_index < number_around:
                chain = chain + pad
            sequence = chain[monomer_index - number_around : monomer_index + number_around + 1]
        return sequence

    def get_position(self, index, frame):

        return frame.particles.position[index]

    def calc_local_density_distribution(self, frame_number, local_size=6):

        frame = self.trajectory.read_frame(frame_number)

        self.get_box(frame_number)
        num_boxes = np.ceil(np.divide(self.box, local_size))
        num_boxes = [int(i) for i in num_boxes]

        a = np.zeros(num_boxes)

        for ind, p in enumerate(frame.particles.position):
            if self.get_type(ind, frame)[-1] != "i":
                p = np.add(np.divide(self.box, 2), p)
                a[int(np.floor(p[0]/local_size))][int(np.floor(p[1]/local_size))][int(np.floor(p[2]/local_size))] += 1
        return a

    def calc_local_charge_density_distribution(self, frame_number, local_size=6, count_divalent=False, count_all=False,
                                               count_ion_only_bins=False):

        frame = self.trajectory.read_frame(frame_number)

        self.get_box(frame_number)
        num_boxes = np.ceil(np.divide(self.box, local_size))
        num_boxes = [int(i) for i in num_boxes]

        #print(self.box, num_boxes)
        #quit()

        a = np.multiply(np.ones(num_boxes), -54.3)
        activated = np.multiply(np.ones(num_boxes), 0)

        for ind, p in enumerate(frame.particles.position):
            if (self.get_type(ind, frame)[-1] != "i" or (count_divalent and self.get_type(ind, frame)[1] == "2") or count_all) and self.get_type(ind, frame)[0] == "Q":
                #print(self.get_type(ind, frame))
                p = np.add(np.divide(self.box, 2), p)
                to_add = 0
                if "M" in self.get_type(ind, frame):
                    to_add = -1
                    if "2" in self.get_type(ind, frame):
                        to_add = -2
                elif "P" in self.get_type(ind, frame):
                    to_add = 1
                    if "2" in self.get_type(ind, frame):
                        to_add = 2
                box = a[int(np.floor(p[0]/local_size))][int(np.floor(p[1]/local_size))][int(np.floor(p[2]/local_size))]
                act = activated[int(np.floor(p[0]/local_size))][int(np.floor(p[1]/local_size))][int(np.floor(p[2]/local_size))]
                if box == -54.3:
                    a[int(np.floor(p[0] / local_size))][int(np.floor(p[1] / local_size))][int(np.floor(p[2] / local_size))] = 0
                a[int(np.floor(p[0]/local_size))][int(np.floor(p[1]/local_size))][int(np.floor(p[2]/local_size))] += to_add
                if act == 0 and self.get_type(ind, frame)[-1] != "i":
                    activated[int(np.floor(p[0]/local_size))][int(np.floor(p[1]/local_size))][int(np.floor(p[2]/local_size))] = 1

        if not count_ion_only_bins:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    for k in range(a.shape[2]):
                        if not activated[i][j][k]:
                            a[i][j][k] = -54.3

        return a


    def calc_local_composition(self, frame_number, local_size):

        frame = self.trajectory.read_frame(frame_number)

        self.get_box(frame_number)
        num_boxes = np.ceil(np.divide(self.box, local_size))
        num_boxes = [int(i) for i in num_boxes]

        # print(self.box, num_boxes)
        # quit()

        a = np.multiply(np.ones(num_boxes), 0)
        h = np.multiply(np.ones(num_boxes), 0)
        tracked_monomers = ["EHMA", "BC"]

        for chain_index in range(len(self.all_chain_indices)):
            backbone_indices_chainx = [self.all_chain_indices[chain_index][i][0] for i in
                                       range(len(self.all_chain_indices[chain_index]))]
            for mon, overall_ind in enumerate(backbone_indices_chainx):
                p = frame.particles.position[overall_ind]
                p = np.add(p, np.divide(self.box, 2))
                #print(p, int(np.floor(p[0] / local_size)))
                a[int(np.floor(p[0] / local_size))][int(np.floor(p[1] / local_size))][
                    int(np.floor(p[2] / local_size))] += 1
                if self.all_chain_sequences[chain_index][mon] in tracked_monomers:
                    h[int(np.floor(p[0] / local_size))][int(np.floor(p[1] / local_size))][
                        int(np.floor(p[2] / local_size))] += 1
        #quit()
        s = a.shape
        for i in range(s[0]):
            for j in range(s[1]):
                for k in range(s[2]):
                    if a[i][j][k] == 0:
                        a[i][j][k] = 1

        return np.divide(h, a), a

    def calc_local_composition2(self, frame_number, local_size, offset_x=0, offset_y=0, offset_z=0):

        frame = self.trajectory.read_frame(frame_number)

        self.get_box(frame_number)
        num_boxes = np.ceil(np.divide(self.box, local_size))
        num_boxes = [int(i) for i in num_boxes]

        # print(self.box, num_boxes)
        # quit()

        a = np.multiply(np.ones(num_boxes), 0)
        h = np.multiply(np.ones(num_boxes), 0)
        tracked_monomers = ["EHMA", "BC"]

        for chain_index in range(len(self.all_chain_indices)):
            backbone_indices_chainx = [self.all_chain_indices[chain_index][i][0] for i in
                                       range(len(self.all_chain_indices[chain_index]))]
            for mon, overall_ind in enumerate(backbone_indices_chainx):
                p = frame.particles.position[overall_ind]
                p = np.add(p, np.divide(self.box, 2))
                p = np.add(p, [offset_x, offset_y, offset_z])
                p = np.remainder(p, self.box)
                a[int(np.floor(p[0] / local_size))][int(np.floor(p[1] / local_size))][
                    int(np.floor(p[2] / local_size))] += 1
                if self.all_chain_sequences[chain_index][mon] in tracked_monomers:
                    h[int(np.floor(p[0] / local_size))][int(np.floor(p[1] / local_size))][
                        int(np.floor(p[2] / local_size))] += 1
        #quit()
        s = a.shape
        for i in range(s[0]):
            for j in range(s[1]):
                for k in range(s[2]):
                    if a[i][j][k] == 0:
                        a[i][j][k] = 1

        return np.divide(h, a), a

    def graph_local_composition(self, frames, local_sizes, save_name="compositions.png", min_monomers=5):

        bes = []
        hists = []
        for local_size in local_sizes:
            all_data = []
            for frame in frames:
                ldda, a = self.calc_local_composition(frame, local_size)
                ldda = ldda[:-1, :-1, :-1]
                a = a[:-1, :-1, :-1]
                ldda = ldda.flatten()
                a = a.flatten()
                keep = [comp for ind, comp in enumerate(ldda) if a[ind] >= min_monomers]
                all_data.extend(keep)


            #all_data = [data for data in all_data if data != -54.3]
            #print(all_data)
            #all_data = np.abs(all_data)
            #bins = np.array(range(int(np.min(all_data) - 1), int(np.max(all_data) + 1), 1))
            #bins = np.subtract(bins, -.0000001)
            hist, be = np.histogram(all_data, bins=np.multiply(.01, list(range(0,101, 5))))
            hist = np.divide(hist, np.sum(hist))

            print(np.sum([hist[i] * (be[i] + be[i + 1]) / 2 for i in range(len(hist))]))
            #be = np.divide(be, np.power(local_size, 3))
            bes.append(be)
            hists.append(hist)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_ylim([0, .5])
        ax1.set_xlim([0, 1])
        ax1.set_xlabel('$F_{H}$')
        ax1.set_ylabel('P($F_{H}$)')
        for i in range(len(hists)):
            hist = hists[i]
            be = bes[i]
            ax1.plot(be[1:], hist, label= "L = " + str(local_sizes[i]* .5) + " nm")
        plt.legend(prop={'size': 13})
        plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_local_composition2(self, frames, local_sizes, save_name="compositions.png", min_monomers=5):

        bes = []
        hists = []
        for local_size in local_sizes:
            all_data = []
            num_offsets = int(np.floor(local_size))
            for frame in frames:
                for i in range(0, num_offsets, int(np.floor(num_offsets/4))):
                    for j in range(0, num_offsets, int(np.floor(num_offsets/4))):
                        for k in range(0, num_offsets, int(np.floor(num_offsets/4))):
                            ldda, a = self.calc_local_composition2(frame, local_size, offset_x=i, offset_y=j, offset_z=k)
                            ldda = ldda[:-1, :-1, :-1]
                            a = a[:-1, :-1, :-1]
                            ldda = ldda.flatten()
                            a = a.flatten()
                            keep = [comp for ind, comp in enumerate(ldda) if a[ind] >= min_monomers]
                            all_data.extend(keep)
                            del a
                            del ldda


            hist, be = np.histogram(all_data, bins=np.multiply(.01, list(range(0,101, 5))))
            hist = np.divide(hist, np.sum(hist))
            print(np.sum([hist[i] * (be[i] + be[i + 1]) / 2 for i in range(len(hist))]))
            bes.append(be)
            hists.append(hist)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_ylim([0, .5])
        ax1.set_xlim([0, 1])
        ax1.set_xlabel('$F_{H}$', size=44)
        ax1.set_ylabel('P($F_{H}$)')
        colors = [(1,0,0), (.7,.3,0), (0, .3, .7), (0, 0, 1)]
        for i in range(len(hists)):
            hist = hists[i]
            be = bes[i]
            ax1.plot(be[1:], hist, label= "L=" + str(local_sizes[i]* .5) + "nm", color=colors[i])
        plt.legend(prop={'size': 13})
        #plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.savefig("compositions_diverging", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_local_density_distribution(self, frames, local_size=6):

        all_data = []
        for frame in frames:
            ldda = self.calc_local_charge_density_distribution(frame, local_size=local_size)
            ldda = ldda[:-1, :-1, :-1]
            ldda = ldda.flatten()
            #print(ldda)
            all_data.extend(ldda)

        bin_number=10
        hist, be = np.histogram(all_data, bins=bin_number)
        be = np.divide(be, np.power(local_size, 3))
        print(be)
        hist = np.divide(hist, np.sum(hist))
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Density')
        ax1.set_ylabel('P(density)')
        ax1.plot(be[1:], hist)
        plt.savefig("local_" + str(local_size)+ "_bins"+ str(bin_number) + ".png", bbox_inches='tight', pad_inches=.2)
        #plt.show()

    def graph_local_sizes(self, frames, local_sizes):

        bes = []
        hists = []
        for local_size in local_sizes:
            all_data = []
            for frame in frames:
                ldda = self.calc_local_charge_density_distribution(frame, local_size=local_size)
                ldda = ldda[:-1, :-1, :-1]
                ldda = ldda.flatten()
                # print(ldda)
                all_data.extend(ldda)

            bin_number = 15
            bins = np.array(range(int(np.min(all_data) - 1), int(np.max(all_data) + 1), 1))
            bins = np.subtract(bins, -.0000001)
            hist, be = np.histogram(all_data, bins=bins)
            hist = np.divide(hist, np.sum(hist))
            #be = np.divide(be, np.power(local_size, 3))
            bes.append(be)
            hists.append(hist)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Density')
        ax1.set_ylabel('P(density)')
        for i in range(len(hists)):
            hist = hists[i]
            print(i, np.sum(hist))
            be = bes[i]
            ax1.plot(be[1:], hist, label=str(local_sizes[i]))
        plt.legend()
        plt.savefig("local_" + str(local_size) + "_bins" + str(bin_number) + ".png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def graph_local_abs_charge(self, frames, local_sizes, count_divalent=False, count_all=False, save_name=None,
                               count_ion_only_bins=False):

        bes = []
        hists = []
        for local_size in local_sizes:
            all_data = []
            for frame in frames:
                ldda = self.calc_local_charge_density_distribution(frame, local_size=local_size,
                                                                   count_divalent=count_divalent, count_all=count_all,
                                                                   count_ion_only_bins=count_ion_only_bins)
                ldda = ldda[:-1, :-1, :-1]
                ldda = ldda.flatten()
                # print(ldda)
                all_data.extend(ldda)

            bin_number = 15
            all_data = [data for data in all_data if data != -54.3]
            #print(all_data)
            print(np.sum(all_data)/len(frames))
            #all_data = np.abs(all_data)
            bins = np.array(range(int(np.min(all_data) - 1), int(np.max(all_data) + 1), 1))
            bins = np.subtract(bins, -.0000001)
            hist, be = np.histogram(all_data, bins=bins)
            hist = np.divide(hist, np.sum(hist))
            #be = np.divide(be, np.power(local_size, 3))
            bes.append(be)
            hists.append(hist)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.set_ylim([0, .5])
        ax1.set_xlim([-25,25])
        ax1.set_xlabel('$z_{eff}$')
        ax1.set_ylabel('P($z_{eff}$)')
        for i in range(len(hists)):
            hist = hists[i]
            be = bes[i]
            ax1.plot(be[1:], hist, label= "L = " + str(local_sizes[i]* .5) + " nm")
        plt.legend(prop={'size': 13})
        if save_name is None:
            plt.savefig("localabs_" + str(local_size) + "_bins" + str(bin_number) + ".png", bbox_inches='tight', pad_inches=.2)
        else:
            plt.savefig(save_name, bbox_inches='tight',
                        pad_inches=.2)
        plt.show()

    def all_polymer_dists(self, frame_number):

        frame = self.trajectory.read_frame(frame_number)

        positions = np.array([pos for ind, pos in enumerate(frame.particles.position)
                              if self.get_type(ind, frame)[:-1] != 'i'])

        print(positions.shape)
        matrix = distances.self_distance_array(positions)
        #matrix = np.triu(matrix)
        dists = matrix.flatten()
        #dists = np.array([data for data in dists if data != 0])
        print("did")
        return dists

    def all_backbone_dists(self, frame_number):

        backbond_indices = [self.all_chain_indices[chain][x][0] for chain in range(len(self.all_chain_indices))
                            for x in range(len(self.all_chain_indices[chain]))]

        frame = self.trajectory.read_frame(frame_number)

        positions = np.array([pos for ind, pos in enumerate(frame.particles.position)
                              if ind in backbond_indices])

        print(positions.shape)
        matrix = distances.self_distance_array(positions)
        # matrix = np.triu(matrix)
        dists = matrix.flatten()
        # dists = np.array([data for data in dists if data != 0])
        print("did")
        return dists


    def all_hydrophobic_dists(self, frame_number, chain_indices=None):

        frame = self.trajectory.read_frame(frame_number)

        if chain_indices is None:
            chain_indices = list(range(len(self.all_chain_indices)))


        positions = np.array([pos for ind, pos in enumerate(frame.particles.position)
                              if self.get_type(ind, frame) == "B"])

        print(positions.shape)
        print(self.all_chain_indices[0])

        quit()
        matrix = distances.self_distance_array(positions)
        # matrix = np.triu(matrix)
        dists = matrix.flatten()
        # dists = np.array([data for data in dists if data != 0])
        print("did")
        return dists

    def n_n_dist(self, frame_number, nncut=1.2):

        frame = self.trajectory.read_frame(frame_number)
        all_hydro_positions = np.array([pos for ind, pos in enumerate(frame.particles.position)
                              if self.get_type(ind, frame) == "B"])

        all_nns = []
        string_indices = []
        surface_indices = []
        glob_indices = []
        for chain_index in range(len(self.all_chain_indices)):
        #for chain_index in [0, 1]:
            backbone_indices_chainx = [self.all_chain_indices[chain_index][i][0] for i in range(len(self.all_chain_indices[chain_index]))]
            backbone_positions = np.array([frame.particles.position[x] for x in backbone_indices_chainx])
            all_distances = distances.distance_array(backbone_positions, all_hydro_positions, box=[55,55,55, 90, 90, 90])
            for i in range(len(backbone_indices_chainx)):
                first = np.array(all_distances[i])
                #for q in first:
                #    print(q)
                #quit()
                count = first[first < nncut].shape[0] - 1
                all_nns.append(count)
                if count < 4:
                    string_indices.append((chain_index, i))
                    print(chain_index, i, "string", self.all_chain_sequences[chain_index][i])
                elif count > 5:
                    glob_indices.append((chain_index, i))
                    print(chain_index, i, "glob", self.all_chain_sequences[chain_index][i])
                else:
                    surface_indices.append((chain_index, i))
                    print(chain_index, i, "surface", self.all_chain_sequences[chain_index][i])


        monomers = ["EHMA", "SPMMA2", "PEGMEMA","BC", "QVP"]
        string_comp = [0,0,0,0,0]
        surface_comp = [0, 0, 0, 0, 0]
        glob_comp = [0, 0, 0, 0, 0]
        for thing in surface_indices:
            x = monomers.index(self.all_chain_sequences[thing[0]][thing[1]])
            surface_comp[x] += 1

        for thing in string_indices:
            x = monomers.index(self.all_chain_sequences[thing[0]][thing[1]])
            string_comp[x] += 1

        for thing in glob_indices:
            x = monomers.index(self.all_chain_sequences[thing[0]][thing[1]])
            glob_comp[x] += 1

        not_string_comp = np.add(surface_comp, glob_comp)

        glob_comp = np.divide(glob_comp, np.sum(glob_comp))
        print("glob", monomers, glob_comp)
        surface_comp = np.divide(surface_comp, np.sum(surface_comp))
        print("surface", monomers, surface_comp)

        print("strings=" + str(np.sum(string_comp)) + ", not_string=" + str(np.sum(not_string_comp)))

        print("without positive  strings=" + str(np.sum(string_comp[:3])) + ", not_string=" + str(np.sum(not_string_comp[:3])))
        string_comp = np.divide(string_comp, np.sum(string_comp))
        print("string", monomers, string_comp)

        not_string_comp = np.divide(not_string_comp, np.sum(not_string_comp))
        print("not_string", monomers, not_string_comp)


        neg_string_comp = string_comp[:3]
        neg_string_comp = np.divide(neg_string_comp, np.sum(neg_string_comp))
        print("neg string", monomers[:3], neg_string_comp)

        neg_not_string_comp = not_string_comp[:3]
        neg_not_string_comp = np.divide(neg_not_string_comp, np.sum(neg_not_string_comp))
        print("neg not string", monomers[:3], neg_not_string_comp)
        quit()

        hist, be = np.histogram(all_nns, bins=np.max(all_nns))
        hist = np.divide(hist, np.sum(hist))
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Nearest Neighbors')
        ax1.set_ylabel('P(N.N.)')
        print(be)
        be2 = list(range(1, len(be)))
        print(be2)
        print(be2[:3], np.sum(hist[:3]))
        #ax1.plot(be[1:], hist)
        ax1.plot(be2, hist)
        plt.savefig("neighbor_dist.png", bbox_inches='tight', pad_inches=.2)
        plt.show()


    def n_n_dist_cum(self, frame_number, nncut=1.2):

        frame = self.trajectory.read_frame(frame_number)
        all_hydro_positions = np.array([pos for ind, pos in enumerate(frame.particles.position)
                              if self.get_type(ind, frame) == "B"])

        all_nns = []
        for chain_index in range(len(self.all_chain_indices)):
        #for chain_index in [0, 1]:
            backbone_indices_chainx = [self.all_chain_indices[chain_index][i][0] for i in range(len(self.all_chain_indices[chain_index]))]
            backbone_positions = np.array([frame.particles.position[x] for x in backbone_indices_chainx])
            all_distances = distances.distance_array(backbone_positions, all_hydro_positions, box=[55,55,55, 90, 90, 90])
            for i in range(len(backbone_indices_chainx)):
                first = np.array(all_distances[i])
                #for q in first:
                #    print(q)
                #quit()
                count = first[first < nncut].shape[0] - 1
                all_nns.append(count)
                print(chain_index, i, count)


        hist, be = np.histogram(all_nns, bins=np.max(all_nns))
        hist = np.divide(hist, np.sum(hist))
        hist_cum = [np.sum(hist[:i]) for i in range(1, len(hist) + 1)]
        print(hist)
        print(hist_cum)
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('Nearest Neighbors')
        ax1.set_ylabel('P(N.N.)')
        print(be)
        be2 = list(range(1, len(be)))
        print(be2)
        print(be2[:3], np.sum(hist[:3]))
        #ax1.plot(be[1:], hist)
        ax1.plot(be2, hist_cum)
        plt.savefig("neighbor_dist_cum.png", bbox_inches='tight', pad_inches=.2)
        plt.show()




    def g_r(self, frames):

        #all_data =self.all_polymer_dists(frames[0])
        all_data = self.all_polymer_dists(frames[0])
        #for frame in frames[1:]:
        #    print(frame)
        #    all_data = np.concatenate(self.all_polymer_dists(frame), all_data)
        #    print("extendeed")


        hist, be = np.histogram(all_data, bins=2000)
        #for i in range(len(hist)):
            #hist[i] = hist[i] / (4/3 * np.pi * np.power(be[i+1], 3))
        #be = np.divide(be, np.power(local_size, 3))

        hist = np.divide(hist, np.sum(hist))
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('R (nm)')
        ax1.set_ylabel('G(r)')
        ax1.plot(be[1:], np.divide(hist, be[1:]**3))
        plt.savefig("gr.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def g_r_return(self, frames):


        all_data = self.all_polymer_dists(frames[0])

        hist, be = np.histogram(all_data, bins=2000)
        del all_data
        hist = np.divide(hist, np.sum(hist))

        bin_middles = [(be[i] + be[i+1])/2 for i in range(len(be) - 1)]

        for ind, frame in enumerate(frames[1:]):
            all_data = self.all_polymer_dists(frame)
            new_hist, new_be = np.histogram(all_data, bins=be)
            del all_data
            new_hist = np.divide(new_hist, np.sum(new_hist))
            hist = np.add(np.multiply(new_hist, 1 / (2 + ind)), np.multiply(hist, (1+ind) / (2+ind)))


        return hist, bin_middles

    def s_q(self, frames, lamdas):

        # all_data =self.all_polymer_dists(frames[0])
        all_data = self.all_backbone_dists(frames[0])
        # for frame in frames[1:]:
        #    print(frame)
        #    all_data = np.concatenate(self.all_polymer_dists(frame), all_data)
        #    print("extendeed")

        qs = np.divide(2 * np.pi, lamdas)

        s_s = []
        for q in qs:
            s_s.append(np.sum([np.sin(q * r) / (q * r) for r in all_data]))
            print(q, "did one")


        fig = plt.figure()

        qs = np.divide(qs, .5)

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('q ($nm^{-1}$)')
        ax1.set_ylabel('S(q)')
        ax1.plot(qs, s_s)
        plt.savefig("sq1.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def s_q2(self, frames, lamdas):

        g, rs = self.g_r_return(frames)
        qs = np.divide(2 * np.pi, lamdas)
        s_s = []
        for q in qs:
            integrand = np.array([np.multiply(g[i], np.sin(q * rs[i]) / (q * rs[i])) for i in range(len(rs))])
            #print("i", integrand)
            #print("t", np.trapz(integrand, x=rs))
            s_s.append(np.sum(integrand))
            print(q, "did one")

        qs = np.divide(qs, .5)
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_xlabel('q ($nm^{-1}$)')
        ax1.set_ylabel('S(q)')
        self.write_file(s_s, qs)
        ax1.plot(qs, s_s)
        plt.savefig("sq22.png", bbox_inches='tight', pad_inches=.2)
        plt.show()

    def write_file(self, s,q, file_name="sq22.txt"):

        f = open(file_name, "w")
        for i in range(len(s)):
            f.write(str(s[i]) + " " + str(q[i]) + "\n")

    def get_type(self, index, frame):

        if not isinstance(index,int):
            #print("Warning index "  + str(index) + "is not an integer")
            return("!!!Faketype")
        return frame.particles.types[frame.particles.typeid[index]]

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

    def condensed_dyes(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed.txt")

        condensed = [dye for dye in self.dyes if self.is_condensed(dye, frame)]

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed.txt", "w")
        for dye in condensed:
            string = str(dye[0])
            string += "\n"
            f.write(string)
        f.close()

        return condensed

    def condensed_dyes_cv(self, frame, always_com=False):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_cv.txt"):
            return self.read_cv(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_cv.txt")

        condensed = [dye for dye in self.dyes if self.is_condensed(dye, frame, cut=1.122)]
        cv = [self.cv(dye, number, always_com=always_com) for dye in condensed]

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_cv.txt", "w")
        for circle, ni in cv:
            string = str(circle) + " " + str(ni)
            string += "\n"
            f.write(string)
        f.close()

        return condensed

    def is_condensed(self, dye, frame, cut=1.122):

        chain_positions = []
        dye_positions = [self.get_position(d, frame) for d in dye if self.get_type(d, frame)[-1] != "i"]

        if len(dye_positions) ==0:
            return False

        for chain in self.rando_chain_indices:
            for mon in chain:
                for bead in mon:
                    if self.get_type(bead, frame)[0] == "B":
                        chain_positions.append(self.get_position(bead, frame))
        #                if len(dye) > 1:
        #                    for d in dye:
        #                        if self.distance(self.get_position(d, frame), self.get_position(bead, frame)) < cut:
        #                            return True
         #               else:
         #                      if self.distance(self.get_position(dye, frame), self.get_position(bead, frame)) < cut:
         #                           return True


        dye_positions = np.array(dye_positions)
        chain_positions = np.array(chain_positions)
        print(dye_positions.shape, chain_positions.shape)
        dist_array = distances.distance_array(dye_positions, chain_positions)
        if np.any(dist_array < cut):
            return True
        return False

    def mon_condenser(self, dye, frame, cut=1.122):

        chain_positions = []
        chain_indices = []
        dye_positions = [self.get_position(d, frame) for d in dye if self.get_type(d, frame)[-1] != "i"]
        if len(dye_positions) ==0:
            return False

        for chain_ind, chain in enumerate(self.all_chain_indices):
            for mon_ind, mon in enumerate(chain):
                for bead in mon:
                    if self.get_type(bead, frame)[0] == "B":
                        chain_positions.append(self.get_position(bead, frame))
                        chain_indices.append((chain_ind, mon_ind))
        #                if len(dye) > 1:
        #                    for d in dye:
        #                        if self.distance(self.get_position(d, frame), self.get_position(bead, frame)) < cut:
        #                            return chain_ind, mon_ind
        #                else:
        #                       if self.distance(self.get_position(dye, frame), self.get_position(bead, frame)) < cut:
        #                            return chain_ind, mon_ind

        dye_positions = np.array(dye_positions)
        chain_positions = np.array(chain_positions)
        print(dye_positions.shape, chain_positions.shape)
        dist_array = distances.distance_array(dye_positions, chain_positions)
        if np.any(dist_array < cut):
            for i in range(dist_array.shape[1]):
                if np.any(dist_array[:, i] < cut):
                    return chain_indices[i]

        return False

    def condensed_dye_sequences(self, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        if os.path.exists(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_seq.txt"):
            return self.read_data(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_seq.txt", seq=True)

        sequences = []

        for dye in self.dyes:
            indices = self.mon_condenser(dye, frame)
            print(indices)
            if indices:
                sequences.append(self.get_sequence_around_monomer(indices[1], indices[0]))

        f = open(str(self.gsd_name[:-4]) + "_frame" + str(number) + "_condensed_seq.txt", "w")
        for sequence in sequences:
            string = ""
            for thing in sequence:
                string += str(thing) + " "
            string += "\n"
            f.write(string)

        f.close()
        return sequences


    def is_condensable(self, dye, frame):

        number = frame
        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        for d in dye:
            if self.get_type(d, frame)[-1] != 'i':
                return True
        return False


    def percent_condensed(self, frames=None):

        data = []

        if frames is None:
            max = len(self.trajectory)
            frames = []
            for i in range(max):
                path = str(self.gsd_name[:-4]) + "_frame" + str(i) + "_condensed_seq.txt"
                if os.path.exists(path):
                    frames.append(i)
                    data.append(len(self.condensed_dye_sequences(i)))
        else:
            for i in frames:
                data.append(len(self.condensed_dye_sequences(i)))

        num_dyes = np.sum([self.is_condensable(dye, frames[0]) for dye in self.dyes])

        total = np.sum(data)

        print("standard deviation: " + str(np.std(data)))
        print("standard error: " + str(np.std(data)/np.sqrt(len(frames) * num_dyes)  ))

        return total / (num_dyes * len(frames))


    def graph_all_condensed(self):

        max = len(self.trajectory)

        frames = []
        data = []

        for i in range(max):
            path = str(self.gsd_name[:-4]) + "_frame" + str(i) + "_condensed_seq.txt"
            if os.path.exists(path):
                frames.append(i)
                num_dyes = np.sum([self.is_condensable(dye, i) for dye in self.dyes])
                data.append(len(self.condensed_dye_sequences(i))/ num_dyes)

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('Frame')
        ax1.set_ylabel('percent condensed')
        ax1.plot(frames, data, label="Average = " + str(self.percent_condensed()))
        plt.legend()
        plt.show()

    def polymer_volume_fraction(self, frame=0):

        types = ['Bb', 'B', 'L', 'QM', 'QP', 'QMsmall']

        diameters = [1.0, 1.0, 1.0, 1.0, 1.0, .6]

        counts = [0, 0, 0, 0, 0, 0]

        self.get_box(frame)
        frame = self.trajectory.read_frame(frame)

        for i in range(len(frame.particles.position)):
            if self.get_type(i, frame) in types:
                counts[types.index(self.get_type(i, frame))] += 1

        vol = 0
        for i in range(5):
            vol += counts[i] * (diameters[i]/2)**3 * 4/3 * np.pi
        frac = vol / (self.box[0] * self.box[1] * self.box[2])
        return frac




    def cv(self, dye, frame_number, cv_cut=6.0, always_com=True):

        frame = self.trajectory.read_frame(frame_number)
        unit_vectors = []
        dye_positions = self.get_positions_chain(dye, frame)
        if len(dye) == 1:
            dye_center = dye_positions[0]
        else:
            dye_positions.pop()
            if always_com:
                dye_center = np.average(dye_positions, axis=0)
            else:
                dye_center = dye_positions.pop()
        connections = self.get_all_connections_from_clusters(frame_number)
        indices = self.mon_condenser(dye, frame, cut=1.122)
        if indices:
            searchable = connections[indices[0]]
            pos2 = self.new_pos(frame.particles.position[self.all_chain_indices[indices[0]][indices[1]][0]], dye_center)
            print(self.distance(dye_center, pos2))
        else:
            searchable = []
        for chain_index in searchable:
            chain = self.all_chain_indices[chain_index]
            positions = [frame.particles.position[i] for mon in chain for i in mon if self.get_type(i, frame) == "B"]
            index2 = [i for mon in chain for i in mon if self.get_type(i, frame) == "B"]
            #print(self.all_chain_indices[indices[0]][indices[1]][0] in index2)
            #print(frame.particles.position[self.all_chain_indices[indices[0]][indices[1]][0]], dye_center)
            #print(self.new_pos(frame.particles.position[self.all_chain_indices[indices[0]][indices[1]][0]], dye_center))
            for pos in positions:
                mon_pos = self.new_pos(pos, dye_center, silence=True)
                if self.distance(dye_center, mon_pos) < cv_cut:
                    unit_vectors.append(np.divide(np.subtract(dye_center, mon_pos), la.norm(np.subtract(dye_center, mon_pos))))
        cv = 1 - la.norm(np.sum(unit_vectors, axis=0))/len(unit_vectors)
        print(cv, len(unit_vectors), searchable, indices)
        return cv, len(unit_vectors)

    def read_cv(self, path):
        f = open(path)
        data = f.readlines()
        out = []
        for line in data:
            s = line.split()
            part = [float(s[0]), int(s[1])]
            out.append(part)
        return out

    def get_all_cvs(self, frames):

        path = str(self.gsd_name[:-4]) + "_frame" + str(frames[0]) + "_condensed_cv.txt"
        cvs = self.read_cv(path)
        for frame in frames[1:]:
            path = str(self.gsd_name[:-4]) + "_frame" + str(frame) + "_condensed_cv.txt"
            cvs.extend(self.read_cv(path))
        return cvs

    def graph_cv_distribution(self, frames):

        the_list = self.get_all_cvs(frames)
        cv = [thing[0] for thing in the_list if thing[1] != 0]
        ni = [thing[1] for thing in the_list if thing[1] != 0]
        #bins = [0, .2, .4, .6, .8, 1]
        bins = [.1 *i for i in range(11)]
        hist, bin_edges = np.histogram(cv, bins=bins)
        hist= np.divide(hist, np.sum(hist))
        print("hist")
        for h in hist:
            print(h)
        print("hist")
        #hist = np.divide(hist, np.sum(hist))

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('Circular Variance')
        ax1.set_ylabel('P(CV)')
        ax1.set_xlim([0,1])
        #ax1.bar(bin_edges[1:], hist, width=.2)
        ax1.plot(bin_edges[1:], hist)
        print(np.sum(hist), len(frames))
        #plt.hist(cv, bins=bins)
        plt.show()

    def cv1(self, frame, chain_index):

        frame_number = frame
        frame = self.trajectory.read_frame(frame)
        chain = self.all_chain_indices[chain_index]
        positions = [frame.particles.position[i] for mon in chain for i in mon if self.get_type(i, frame) == "B"]
        positions = [self.new_pos(pos, positions[0]) for pos in positions]
        cvs = []

        test_cvs = []

        for pos in positions:
            unit_vectors = [np.subtract(pos, thing) for thing in positions if self.distance(pos, thing) != 0
                            and self.distance(pos, thing) < 6]

            unit_vectors = [np.divide(unit_vector, la.norm(unit_vector)) for unit_vector in unit_vectors]
            cv = la.norm(np.divide(np.sum(unit_vectors, axis=0), len(unit_vectors)))
            cvs.append(1 - cv)
        for i in range(len(positions)):
            print("[ ",positions[i],"," ,cvs[i], "]")

        average_point = np.average(positions, axis=0)
        print(average_point)

        fig = plt.figure()

        ax1 = fig.add_subplot(111, projection='3d')

        #ax1.set_title("Frame " + str(frame))
        positions = np.array(positions)
        colors = []
        for cv in cvs:
            if cv < 0.2:
                colors.append('r')
            elif cv < 0.4:
                colors.append('y')
            elif cv < 0.6:
                colors.append('g')
            elif cv < 0.8:
                colors.append('b')
            elif cv < 1.0:
                colors.append('m')

        ax1.scatter(positions[:,0], positions[:,1], zs=positions[:,2], s=100, c=colors)
        #ax1.contour(positions[:,0], positions[:,1], positions[:,2])
        plt.legend()
        plt.show()
        name = "file_" + str(frame_number) + ".xyz"
        f = open(name, 'w')
        f.write(str(len(positions)))
        f.write("\n\n")
        for i in range(len(positions)):
            f.write(str(cvs[i]) + " " + str(positions[i][0]) + " " + str(positions[i][1]) + " " + str(positions[i][2])+ "\n")

    def electrostatic_potential(self, charge_index, frame, lb=1.4):

        #frame = self.trajectory.read_frame(frame_number)


        charged_indices = []

        for ind in range(len(frame.particles.position)):
            if "Q" in self.get_type(ind, frame) and ind != charge_index:
                charged_indices.append(ind)
        #charged_indices = np.array([ind for ind in range(len(frame.particles.position)) if "Q" in self.get_type(ind, frame)] and ind != charge_index)

        charges = []

        for index in charged_indices:
            charge = 1
            if "M" in self.get_type(index, frame):
                charge = -1 * charge
            if "2" in self.get_type(index, frame):
                charge = charge * 2
            charges.append(charge)


        ref = frame.particles.position[charge_index]

        ref_charge = 1
        if "M" in self.get_type(charge_index, frame):
            ref_charge = -1 * ref_charge
        if "2" in self.get_type(charge_index, frame):
            ref_charge = ref_charge * 2

        potential = 0

        for o_ind, ind in enumerate(charged_indices):
            new_pos = self.new_pos(frame.particles.position[ind], ref, silence=True)
            dist = self.distance(new_pos, ref)
            potential += (ref_charge * charges[o_ind] * lb / dist)

        return potential


    def condensed_charged_dye_potentials(self, frame, condensed=True):

        #frame = self.trajectory.read_frame(frame_number)
        chared_indices = []
        eps = []

        for dye in self.dyes:
            dye_types = [self.get_type(ind, frame) for ind in dye]
            if "B" in dye_types and ("QPi" in dye_types or "QMi" in dye_types):
                if self.is_condensed(dye, frame) == condensed:
                    chared_indices.append(dye[-2])

        for ci in chared_indices:
            eps.append(self.electrostatic_potential(ci, frame))

        return eps

    def graph_condensed_dye_ep(self, frames, condensed=True):


        all_eps = []
        averages = []

        for frame_number in frames:
            frame = self.trajectory.read_frame(frame_number)
            eps = self.condensed_charged_dye_potentials(frame, condensed=condensed)
            all_eps.extend(eps)
            averages.append(np.average(eps))

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames[-1]))

        ax1.set_xlabel('Frame')
        ax1.set_ylabel('Average EP')
        #ax1.set_xlim([0,1])
        #ax1.bar(bin_edges[1:], hist, width=.2)
        ax1.plot(frames, averages)
        print(len(all_eps), len(frames))
        print(np.average(all_eps))
        print(np.std(all_eps), np.std(all_eps)/ np.sqrt(len(all_eps)))
        #print(np.sum(hist), len(frames))
        #plt.hist(cv, bins=bins)
        plt.show()
