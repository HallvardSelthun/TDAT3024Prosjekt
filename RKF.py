import numpy as np

class RKFmodel():
    def __init__(self, F, w_0, h, T, t_n):
        self.F = F
        self.w_0 = w_0
        self.h = h
        self.T = T
        self.t_n = t_n

    def s1(self, w_i):
        return self.F(w_i)

    def s2(self, w_i):
        return self.F(w_i + (self.h/4)*self.s1(w_i))

    def s3(self, w_i):
        return self.F(w_i + (3*self.h/32)*self.s1(w_i) + (9*self.h/32)*self.s2(w_i))

    def s4(self, w_i):
        return self.F(w_i + (1932*self.h/2197)*self.s1(w_i) - (7200*h/2197)*self.s2(w_i) + (7296*self.h/2197)*self.s3(w_i))

    def s5(self, w_i):
        return self.F(w_i + (439*self.h/216)*self.s1(w_i) - 8*self.h*self.s2(w_i) + (3680*self.h/513)*self.s3(w_i) - (845*self.h/4104)*self.s4(w_i))

    def s6(self, w_i):
        return self.F(w_i - (8*self.h/27)*self.s1(w_i) + 2*self.h*self.s2(w_i) - (3544*self.h/2565)*self.s3(w_i) + (1859*self.h/4104)*self.s4(w_i) - (11*self.h/40)*self.s5(w_i))

    def w_i_pluss_1(self, w_i):
        return w_i + self.h*(25*self.s1(w_i)/216 + 1408*self.s3(w_i)/2565 + 2197*self.s4(w_i)/4104 - self.s5(w_i)/5)

    def z_i_pluss_1(self, w_i):
        return w_i + self.h*(16*self.s1(w_i)/135 + 6656*self.s3(w_i)/12825 + 28561*self.s4(w_i)/56430 - 9*self.s5(w_i)/50 + 2*self.s6(w_i)/55)

    def e_i_pluss_1(self, w_i):
        return np.abs(self.w_i_pluss_1(w_i) + self.z_i_pluss_1(w_i))

    def w_n(self):
        w_i = self.w_0
        t_i = 0

        while t_i<=self.t_n:
            w_i = self.w_i_pluss_1(w_i)
            t_i += self.h

        return w_i

rfk = RKFmodel()