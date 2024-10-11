from DyeAbs import DyeAbs
import numpy as np

class NeutralViolet(DyeAbs):

    def __init__(self):

        super(NeutralViolet, self).__init__(index=-1)

    def build(self, with_ion=False):

        self.position.append([0, 0, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.position.append([-1, 0, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.position.append([0.5, np.sqrt(3)/2, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.position.append([0.5, -np.sqrt(3)/2, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.bonds.append([0, 1])
        self.bond_names.append("polybond")

        self.bonds.append([0, 2])
        self.bond_names.append("polybond")

        self.bonds.append([0, 3])
        self.bond_names.append("polybond")

        self.angles.append([1, 0, 2])
        self.angle_names.append("stiff_3")

        self.angles.append([3, 0, 1])
        self.angle_names.append("stiff_3")

        self.angles.append([2, 0, 3])
        self.angle_names.append("stiff_3")