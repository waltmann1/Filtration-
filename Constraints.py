from __future__ import division
import numpy as np


class Constraints(object):

    def __init__(self):

        self.names = []
        self.r0 = []

        self.names.append("STY-STY")
        self.r0.append(.27)

        self.names.append("SCY-STY")
        self.r0.append(.22)

        self.funct = [1 for _ in range(len(self.names))]

    def martini_string(self, constraint, name):

        constraint = np.add(constraint, 1)
        b0 = self.r0[self.names.index(name)]
        funct = self.funct[self.names.index(name)]
        arr = [constraint[0], constraint[1], funct, b0]
        tring = ""
        for i in range(3):
            tring += str(arr[i]) + " "
        return tring + str(b0) + "\n"
