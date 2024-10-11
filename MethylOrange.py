from DyeAbs import DyeAbs


class MethylOrange(DyeAbs):

    def __init__(self):

        super(MethylOrange, self).__init__(index=-1)

    def build(self, with_ion=True):

        self.position.append([0, 0, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.position.append([1, 0, 0])
        self.charge.append(0)
        self.mass.append(1)
        self.type.append('B')

        self.position.append([-1, 0, 0])
        self.charge.append(-1)
        self.mass.append(1)
        self.type.append('QM')

        if with_ion:
            self.position.append([-2, 0, 0])
            self.charge.append(1)
            self.mass.append(1)
            self.type.append('QPi')

        self.bonds.append([0, 1])
        self.bond_names.append("polybond")

        self.bonds.append([0, 2])
        self.bond_names.append("polybond")

        self.angles.append([1, 0, 2])
        self.angle_names.append("stiff")