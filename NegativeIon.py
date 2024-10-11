from DyeAbs import DyeAbs


class NegativeIon(DyeAbs):

    def __init__(self):

        super(NegativeIon, self).__init__(index=-1)

    def build(self, with_ion=True):

        self.position.append([0, 0, 0])
        self.charge.append(-1)
        self.mass.append(1)
        self.type.append('QMi')

