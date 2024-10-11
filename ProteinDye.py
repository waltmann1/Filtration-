from DyeAbs import DyeAbs
import numpy as np


class ProteinDye(DyeAbs):

    def __init__(self, file_name):

        self.file_name = file_name
        super(ProteinDye, self).__init__(index=0)

    def build(self, with_ion=True):

        f = open(self.file_name)
        data = f.readlines()
        data = data[2:]
        for line in data:
            s = line.split()
            self.position.append([float(s[1]), float(s[2]), float(s[3])])
            self.type.append(s[0])
            if s[0] == "QM":
                self.charge.append(-1)
            elif s[0] == "QP":
                self.charge.append(1)
            else:
                self.charge.append(0)
            self.mass.append(1)

        self.position = np.array(self.position)
        self.position = np.subtract(self.position, np.average(self.position, axis=0))
        self.position = list(self.position)
        # for ind,pos in enumerate(self.position):
        #    for ind2, pos2 in enumerate(self.position[ind + 1:]):
        #       dist = np.linalg.norm(np.subtract(pos, pos2))
        #      if dist < 15:
        #         self.bond_names.append("GNM_bond_" + str(round(dist))[0])
        #        self.bonds.append([ind, ind + 1 + ind2])

        radius = np.max([np.linalg.norm(pos) for pos in self.position])
        if with_ion:
            for i in range(len(self.type)):
                if self.type[i] == "QM":
                    self.type.append("QPi")
                    self.charge.append(1)
                    self.position.append(np.multiply(self.position[i], (radius + 2) / np.linalg.norm(self.position[i])))
                    self.mass.append(1)
                elif self.type[i] == "QP":
                    self.type.append("QMi")
                    self.charge.append(-1)
                    self.mass.append(1)
                    self.position.append(np.multiply(self.position[i], (radius + 2 / np.linalg.norm(self.position[i]))))

        self.position = np.divide(self.position, 5)
