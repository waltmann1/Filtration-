from MonomerAbs import MonomerAbs
import numpy as np


class BC(MonomerAbs):

    def __init__(self, with_ion=False):
        super(BC, self).__init__()
        self.length = 0
        self.mass = 0

    def add_to_polymer(self, p, spot, up=True):

        p.charge[spot] = 0
        p.type[spot] = 'B'