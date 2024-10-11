from MonomerAbs import MonomerAbs
import numpy as np
from MartiniPS import MartiniPS


class Martini2qVP(MartiniPS):

    def __init__(self, with_ion=False):

        super(Martini2qVP, self).__init__()

    def add_to_polymer(self, p, spot, up=True):

        super(Martini2qVP, self).add_to_polymer(p, spot, up=up)
        p.type[spot] = "P4"
        p.type[p.num_beads - 3] = "Qd"
        p.charge[p.num_beads - 3] = 1