from __future__ import division
import numpy as np


class PositionRestraints(object):

    def __init__(self):

        self.names = []
        self.ks = []

        self.names.append("normal")
        self.ks.append([1000, 1000, 1000])

        self.names.append("x_y_only")
        self.ks.append([1000, 1000, 0])

        self.names.append("times10")
        self.ks.append([10000, 10000, 10000])

        self.funct = [1 for _ in range(len(self.names))]

    def martini_string(self, restraint, name):

        restraint = np.add(restraint, 1)
        ks = self.ks[self.names.index(name)]
        funct = self.funct[self.names.index(name)]
        arr = [restraint, funct, ks[0], ks[1], ks[2]]
        tring = ""
        for i in range(4):
            tring += str(arr[i]) + " "
        return tring + str(ks[2]) + "\n"
