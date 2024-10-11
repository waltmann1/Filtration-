from __future__ import division
import numpy as np

from Loggable import Loggable


class Bonds(Loggable):

    def __init__(self, log_list=None):

        super(Bonds, self).__init__(log_list)
        #self.log_values = ['bond_fene_energy']
        self.log_values = ['bond_harmonic_energy']
        self.names = []
        self.k = []
        self.r0 = []
        self.bond_ref = None

        #self.k.append(12)
        self.k.append(120)
        #self.r0.append(1.0)
        self.r0.append(1.0)
        self.names.append('polybond')

        #self.k.append(150000)
        self.k.append(1500)
        self.r0.append(.38)
        self.names.append("go_bond")

        self.k.append(1500)
        self.r0.append(.47)
        self.names.append("flank")


        #self.k.append(12)
        self.k.append(120)
        self.r0.append(1.0)
        self.names.append('sidechain')
        for i in range(16):
            self.r0.append(i * .1)
            self.names.append("GNM_bond_" + str(i))
            self.k.append(10)
        #self.k.append(17000)
        #self.r0.append(0.5)

        self.names.append('STY-STY')
        self.k.append(17000)
        self.r0.append(0.27)

        self.names.append('STY-Na_1')
        self.k.append(17000)
        self.r0.append("PET")

        self.names.append('STY-Na_2')
        self.k.append(17000)
        self.r0.append("PET")

        self.names.append('Na-STY')
        self.k.append(17000)
        self.r0.append("PET")

        self.names.append("Na-Na")
        self.k.append(17000)
        self.r0.append("PET")

        self.names.append("MMABackbone")
        self.k.append(21100)
        self.r0.append(.289)

        self.names.append("SC1-Na")
        self.k.append(17000)
        self.r0.append(.282)

        self.names.append("peg")
        self.k.append(17000)
        self.r0.append(.33)

        self.names.append("peo")
        self.k.append(7000)
        self.r0.append(.32)

        self.names.append("exclusion")
        self.k.append(0)
        self.r0.append(0)

        self.names.append("EHMA1")
        self.k.append(17000)
        self.r0.append(0.54)

        self.names.append("EHMA2")
        self.k.append(1250)
        self.r0.append(0.425)

        self.names.append("SCY-SCY")
        self.k.append(8000)
        self.r0.append(.25)

        self.funct = [6 for _ in range(len(self.k))]


    def set_all_bonds(self, system, attach=0.38, k_att=1500, all_harm=False):
        """

        :param system: the system that needs the parameters set
        :return: reference to the harmonic bond object
        """

        import hoomd
        from hoomd import md
        #self.bond_ref = hoomd.md.bond.fene()

        snap = system.take_snapshot(all=True)
        first = snap.bonds.types[0]
        switched = False
        already_set = []
        if first[:2] == "go":
            self.bond_ref = hoomd.md.bond.table(width=10000)
            self.log_values = ["bond_table_energy"]
            switched = True
        else:
            self.bond_ref = hoomd.md.bond.harmonic()
            self.log_values = ["bond_harmonic_energy"]

        for b in snap.bonds.types:
            name = str(b)
            print(b, b[:2])
            if switched:
                s = name.split("_")
                print(s)
                #quit()
                if s[1] == "morse" and name not in already_set:
                    r0 = float(s[2])
                    d = float(s[3])
                    beta = float(s[4])
                    cut = 6 * r0
                    cut=10
                    func = self.gro_morse
                    if not all_harm:
                        self.bond_ref.bond_coeff.set(name, func=func, rmin=0.1, rmax=cut, coeff=dict(d0=d, beta=beta, r0=r0))
                    else:
                        func= self.harmonic
                        self.bond_ref.bond_coeff.set(name, func=func, rmin=0.1, rmax=system.box.Lx, coeff=dict(r0=r0, k=1500))
                    print(d, beta, r0)
                    already_set.append(name)
                elif s[1] == "harm" and name not in already_set:
                    r0 = float(s[2])
                    k = float(s[3])
                    #print(r0, k)
                    #quit()
                    cut = 6 * r0
                    cut=10
                    func = self.harmonic
                    self.bond_ref.bond_coeff.set(name, func=func, rmin=0.1, rmax=system.box.Lx, coeff=dict(r0=r0, k=k))
                    already_set.append(name)
                elif s[1] == "attach" and name not in already_set:
                    func = self.harmonic
                    r0 = attach
                    k = k_att
                    cut = 20
                    self.bond_ref.bond_coeff.set(name, func=func, rmin=0.1, rmax=cut, coeff=dict(r0=r0, k=k))
                    already_set.append(name)
                else:
                    if s[1] == "attach":
                        print("miss")
                        quit()
                    sigma = float(s[2])
                    epsilon = float(s[3])
                    cut = 3 * sigma
                    cut=10
                    func = self.lj
                    self.bond_ref.bond_coeff.set(name, func=func, rmin=0.1, rmax=cut, coeff=dict(epsilon=epsilon, sigma=sigma))
                    already_set.append(name)
            #self.bond_ref.bond_coeff.set(name, k=self.k[self.names.index(name)], r0=self.r0[self.names.index(name)]
                                        # , epsilon=0.5, sigma=0.3)
            if not switched:
                self.bond_ref.bond_coeff.set(name, k=self.k[self.names.index(name)], r0=self.r0[self.names.index(name)])
        self.add_to_logger()

        return self.bond_ref

    def read_go_name(self, name):

        s = name.split()

    def martini_string(self, bond, name):

        bond = np.add(bond, 1)

        if "PET" in self.r0:
            from MartiniPET import MartiniPET
            mp = MartiniPET()
            m_r0s, names = mp.bonds()
            print(m_r0s, names)
            for ind, theta in enumerate(self.r0):
                if theta == "PET":
                    self.r0[ind] = m_r0s[names.index(self.names[ind])]

        k = self.k[self.names.index(name)]
        b0 = self.r0[self.names.index(name)]
        funct = self.funct[self.names.index(name)]

        arr = [bond[0], bond[1], funct, b0, k]
        tring = ""
        for i in range(4):
            tring += str(arr[i]) + " "

        return tring + str(k) + "; " + name + "\n"

    def gro_morse(self, r, rmin, rmax, d0, beta, r0):

        #print(rmin, rmax, d0, beta, r0)

        V = d0 * ((1 - np.exp(-beta * (r - r0)))**2)
        F = - 2 * d0 * beta * np.exp(- beta * (r - r0)) * (1 - np.exp(-beta * (r - r0)))
        return (V, F)

    def lj(self, r, rmin, rmax, epsilon, sigma):

        V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
        F = 4 * epsilon / r * (12 * (sigma / r)**12 - 6 * (sigma / r)**6)

        return (V, F)

    def harmonic(self, r, rmin, rmax, k, r0):

        V = .5 * k * (r - r0)**2
        F = - k * (r - r0)

        return (V, F)






