from __future__ import division
import numpy as np

from Loggable import Loggable


class Dihedrals(Loggable):

    def __init__(self, log_list=None):

        super(Dihedrals, self).__init__(log_list)
        #self.log_values = ['bond_fene_energy']
        self.log_values = ['dihedral_harmonic_energy']
        self.names = []
        self.k = []
        self.theta0 = []
        self.dih_ref = None

        self.k.append(50)
        self.theta0.append(0)
        self.names.append('improper')

        self.k.append(167)
        self.theta0.append(180)
        self.names.append('Felipe')


        self.funct = [2 for _ in range(len(self.names))] # improper dihedral

    def set_all_dihedrals(self):
        import hoomd
        from hoomd import md


    def martini_string(self, dih, name):

        dih = np.add(dih, 1)
        k = self.k[self.names.index(name)]
        t0 = self.theta0[self.names.index(name)]
        funct = self.funct[self.names.index(name)]

        arr = [dih[0], dih[1], dih[2], dih[3], funct, t0, k]
        tring = ""
        for i in range(6):
            tring += str(arr[i]) + " "
        return tring + str(k) + "\n"