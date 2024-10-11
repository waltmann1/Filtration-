from DyeAbs import DyeAbs


class PositiveIon(DyeAbs):

    def __init__(self):

        super(PositiveIon, self).__init__(index=-1)

    def build(self, with_ion=False):


        self.position.append([0, 0, 0])
        self.charge.append(1)
        self.mass.append(1)
        self.type.append('QPi')
