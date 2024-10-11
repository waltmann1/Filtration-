from DyeAbs import DyeAbs


class DoubleDivalentIons(DyeAbs):

    def __init__(self):

        super(DoubleDivalentIons, self).__init__(index=-1)

    def build(self, with_ion=True):

        self.position.append([0, 0, 0])
        self.charge.append(-1)
        self.mass.append(1)
        self.type.append('Q2Mi')

        self.position.append([-1, 0, 0])
        self.charge.append(2)
        self.mass.append(1)
        self.type.append('Q2Pi')