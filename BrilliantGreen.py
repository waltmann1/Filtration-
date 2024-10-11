from DyeAbs import DyeAbs
import numpy as np

class BrilliantGreen(DyeAbs):

    def __init__(self):

        super(BrilliantGreen, self).__init__(index=-1)

    def build(self, with_ion=True):

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

        self.position.append([1.0, np.sqrt(3), 0])
        self.charge.append(0.5)
        self.mass.append(1)
        self.type.append('QP')

        self.position.append([1.0, -np.sqrt(3), 0])
        self.charge.append(0.5)
        self.mass.append(1)
        self.type.append('QP')

        if with_ion:
            self.position.append([1/2, np.sqrt(3) + 1, 0])
            self.charge.append(-1)
            self.mass.append(1)
            self.type.append('QMi')

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

        self.bonds.append([2, 4])
        self.bond_names.append("polybond")

        self.bonds.append([3, 5])
        self.bond_names.append("polybond")

        self.angles.append([0, 2, 4])
        self.angle_names.append("stiff")

        self.angles.append([0, 3, 5])
        self.angle_names.append("stiff")