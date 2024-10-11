from DyeAbs import DyeAbs


class DivalentNegativeIon(DyeAbs):

    def __init__(self):

        super(DivalentNegativeIon, self).__init__(index=-1)

    def build(self, with_ion=False):

        self.position.append([0, 0, 0])
        self.charge.append(-2)
        self.mass.append(1)
        self.type.append('Q2Mi')
