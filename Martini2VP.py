from MonomerAbs import MonomerAbs
import numpy as np
from MartiniPS import MartiniPS


class Martini2VP(MartiniPS):

    def __init__(self, with_ion=False):

        super(Martini2VP, self).__init__()

    def add_to_polymer(self, p, spot, up=True):

        super(Martini2VP, self).add_to_polymer(p, spot, up=up)
        p.type[spot] = "P4"
        p.type[p.num_beads - 3] = "P4"
