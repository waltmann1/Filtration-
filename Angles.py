from __future__ import division
import numpy as np

from Loggable import Loggable


class Angles(Loggable):

    def __init__(self, log_list=None):

        super(Angles, self).__init__(log_list)
        self.log_values = ['angle_harmonic_energy']
        self.names = []
        self.k = []
        self.theta = []
        self.angle_ref = None

        #self.names.append('nothing')
        #self.theta.append(np.deg2rad(180))
        #self.k.append(0)

        self.names.append('stiff')
        self.theta.append(np.deg2rad(180))
        self.k.append(10000)

        self.names.append('stiff90')
        self.theta.append(np.deg2rad(90))
        self.k.append(10000)

        self.names.append('stiff_3')
        self.theta.append(np.deg2rad(120))
        self.k.append(10000)

        self.names.append('k1')
        self.theta.append(np.deg2rad(180))
        self.k.append(1)

        self.names.append('k10')
        self.theta.append(np.deg2rad(180))
        self.k.append(10)

        self.names.append('k100')
        self.theta.append(np.deg2rad(180))
        self.k.append(100)

        self.names.append('k5')
        self.theta.append(np.deg2rad(180))
        self.k.append(5)

        self.names.append('Na-Na-STY')
        self.theta.append("PET")
        self.k.append(35)

        self.names.append('Na-STY-STY')
        self.theta.append("PET")
        self.k.append(35)

        self.names.append('STY-STY-STY')
        self.theta.append("PET")
        self.k.append(35)

        self.names.append('STY-STY-Na')
        self.theta.append("PET")
        self.k.append(35)

        self.names.append('STY-Na-Na')
        self.theta.append("PET")
        self.k.append(35)

        self.names.append("MMABackbone")
        self.theta.append(175)
        self.k.append(13)

        self.names.append("mma_join")
        self.theta.append(144)
        self.k.append(67)

        self.names.append("peg")
        self.theta.append(130)
        self.k.append(50)

        self.names.append("peo")
        self.theta.append(122)
        self.k.append(400)

        self.names.append("Na-flex-Na")
        self.theta.append(130)
        self.k.append(50)

        self.names.append("EHMA")
        self.theta.append(180)
        self.k.append(25)

        self.names.append("B-R-R")
        self.theta.append(156)
        self.k.append(200)

        self.names.append("B-B-B")
        self.theta.append(170)
        self.k.append(25)

        self.names.append("R-B-B")
        self.theta.append(125)
        self.k.append(45)

        self.funct = [1 for _ in range(len(self.names))]

    def set_all_harmonic_angles(self, system, reset=False, poly=0, go_angles=None):

        import hoomd

        if reset:
            self.angle_ref.disable()
        if len(system.angles) == 0:
            return

        self.angle_ref = hoomd.md.angle.harmonic()
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for a in snap.angles.types:
            name = str(a)
            if name[:2] == 'go':
                s = name.split('_')
                self.names.append(name)
                if go_angles is not None:
                    self.k.append(go_angles)
                else:
                    self.k.append(float(s[3]))
                self.theta.append(np.deg2rad(float(s[2])))
            self.angle_ref.angle_coeff.set(name, k=self.k[self.names.index(name)],
                                           t0=self.theta[self.names.index(name)])
        return self.angle_ref

    def martini_string(self, angle, name):

        angle = np.add(angle, 1)
        if "PET" in self.theta:
            from MartiniPET import MartiniPET
            mp = MartiniPET()
            m_thetas, names = mp.angles()
            print(m_thetas, names)
            for ind, theta in enumerate(self.theta):
                if theta == "PET":
                    self.theta[ind] = m_thetas[names.index(self.names[ind])]

        k = self.k[self.names.index(name)]
        theta = self.theta[self.names.index(name)]
        funct = self.funct[self.names.index(name)]

        arr = [angle[0], angle[1], angle[2], funct, theta, k]
        tring = ""
        for i in range(5):
            tring += str(arr[i]) + " "
        return tring + str(k) + "\n"



