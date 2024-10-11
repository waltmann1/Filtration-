from DyeAbs import DyeAbs


class OneBeadDye(DyeAbs):

    def __init__(self):

        super(OneBeadDye, self).__init__(index=-1)

    def build(self):

        self.position.append([0, 0, 0])
        self.charge.append(0)
        self.mass.append(3)
        self.type.append('B')
