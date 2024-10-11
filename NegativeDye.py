from DyeAbs import DyeAbs


class NegativeDye(DyeAbs):

    def __init__(self):

        super(NegativeDye, self).__init__(index=0)

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