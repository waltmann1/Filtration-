from PEGMEMA import PEGMEMA
import numpy as np


class OneMEMA(PEGMEMA):

    def __init__(self):
        super(OneMEMA, self).__init__(n=1)
