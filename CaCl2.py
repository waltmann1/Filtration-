from DyeAbs import DyeAbs


class CaCl2(DyeAbs):

    def __init__(self):

        super(CaCl2, self).__init__(index=-1)

    def build(self, with_ion=True):

        self.position.append([0, 0, 0])
        self.charge.append(-1)
        self.mass.append(1)
        self.type.append('QMi')

        self.position.append([1, 0, 0])
        self.charge.append(-1)
        self.mass.append(1)
        self.type.append('QMi')

        self.position.append([-1, 0, 0])
        self.charge.append(2)
        self.mass.append(1)
        self.type.append('Q2Pi')
