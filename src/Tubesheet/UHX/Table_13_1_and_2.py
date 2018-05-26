import decimal as d
import math


class _Table_13_1_and_2_base():
    def __init__(self, xa, v_star, q_3=None):
        # See bottom of class for properties
        self.xa = d.Decimal(xa)
        self.v_star = d.Decimal(v_star)
        if q_3 is not None:
            self.__q_3 = d.Decimal(q_3)
        else:
            self.__q_3 = q_3
        self.m = int((4 + (xa / 2) / 2)) * 2 + 1
        d.getcontext().prec = 28  # Python 3.5 default; just enforcing consistency

    # If the bessel functions are too slow it's possible to multi thread with map() and multiprocessing.pool()
    def _ber(self, x):  # Bessel function of the first kind, order 0; Kelvin function
        result = d.Decimal(0)
        if (x <= 0) or (x > self.xa):
            raise ValueError("%s not in range 0 < x < Xa = %s" % (str(x), self.xa))
        if self.m <= 0:
            raise ValueError("m must be greater than 0.")
        for i in range(self.m, -1, -1):
            result = result + ((d.Decimal((-1) ** i)) * (
                (x / d.Decimal(2)) ** (d.Decimal(4.0 * i)) / d.Decimal((math.factorial(2.0 * i)) ** 2.0)))
        return result

    def _bei(self, x):
        result = d.Decimal(0)
        if (x <= 0) or (x > self.xa):
            raise ValueError("%s not in range 0 < x < Xa = %s" % (str(x), self.xa))
        if self.m <= 0:
            raise ValueError("m must be greater than 0.")
        for i in range(self.m, -1, -1):
            result = result + d.Decimal(((-1) ** (i - 1)) * (x / d.Decimal(2)) ** (d.Decimal(4 * i - 2)) / d.Decimal(
                (math.factorial(2 * i - 1)) ** 2))
            if i == 1:
                return result

    def _ber_prime(self, x):
        result = d.Decimal(0)
        if (x <= 0) or (x > self.xa):
            raise ValueError("%s not in range 0 < x < Xa = %s" % (str(x), self.xa))
        if self.m <= 0:
            raise ValueError("m must be greater than 0.")
        for i in range(self.m - 1, -1, -1):
            result += (d.Decimal((-1) ** i)) * d.Decimal(2 * i) * ((x / d.Decimal(2)) ** (d.Decimal((4 * i)) - 1) /
                                                                   d.Decimal((math.factorial(2 * i)) ** 2))
            if i == 1:
                return result

    def _bei_prime(self, x):
        result = d.Decimal(0)
        if (x <= 0) or (x > self.xa):
            raise ValueError("%s not in range 0 < x < Xa = %s" % (str(x), self.xa))
        if self.m <= 0:
            raise ValueError("m must be greater than 0.")
        for i in range(self.m, -1, -1):
            result = result + d.Decimal(((-1) ** (i - 1)) * d.Decimal(2 * i - 1) *
                                        (x / d.Decimal(2)) ** (d.Decimal(4 * i - 3)) /
                                        d.Decimal((math.factorial(2 * i - 1)) ** 2))
            if i == 1:
                return result

    def _phi_sub_one(self, x):
        return self._bei(x) + (d.Decimal((1 - self.v_star)) / x) * self._ber_prime(x)

    def _phi_sub_two(self, x):
        return self._ber(x) - (d.Decimal(1 - self.v_star) / x) * self._bei_prime(x)

    def _calc_Z_a(self):
        return self._bei_prime(self.xa) * self._phi_sub_two(self.xa) - \
               self._ber_prime(self.xa) * self._phi_sub_one(
                   self.xa)  # \ is NOT a division symbol,it's line continuation

    def _calc_Z_d(self):
        num = self._ber(self.xa) * self._phi_sub_two(self.xa) + self._bei(self.xa) * self._phi_sub_one(self.xa)
        den = ((self.xa ** 3) * self._calc_Z_a())
        return num / den

    def _calc_Z_v(self):
        num = self._ber_prime(self.xa) * self._phi_sub_two(self.xa) + \
              self._bei_prime(self.xa) * self._phi_sub_one(self.xa)
        den = (self.xa ** 2) * self._calc_Z_a()
        return num / den

    def _calc_Z_w(self):
        num = self._ber_prime(self.xa) * self._ber(self.xa) + self._bei_prime(self.xa) * self._bei(self.xa)
        den = (self.xa ** 2) * self._calc_Z_a()
        return num / den

    def _calc_Z_m(self):
        num = (self._ber_prime(self.xa) ** 2) + (self._bei_prime(self.xa) ** 2)
        den = self.xa * self._calc_Z_a()
        return num / den

    @property
    def q_3(self):
        return self.__q_3

    @q_3.setter
    def q_3(self, newQ3):
        self.q_3 = d.Decimal(newQ3)

    Z_a = property(fget=_calc_Z_a)
    Z_d = property(fget=_calc_Z_d)
    Z_v = property(fget=_calc_Z_v)
    Z_w = property(fget=_calc_Z_w)
    Z_m = property(fget=_calc_Z_m)


class Table_13_1(_Table_13_1_and_2_base):
    def _calc_Q_m(self, x):
        return (self._bei_prime(self.xa) * self._phi_sub_two(x) - self._ber_prime(self.xa) * self._phi_sub_one(
            x)) / self._calc_Z_a()

    def _calc_Q_v(self, x):
        return (self._phi_sub_one(self.xa) * self._phi_sub_two(x) - self._phi_sub_two(self.xa) * self._phi_sub_one(
            x)) / (self.xa * self._calc_Z_a())

    def _calc_F_m(self, x):
        if self.q_3 is None:
            raise ValueError("q_3 is not set in this instance!")
        else:
            return (self._calc_Q_v(x) + self.q_3 * self._calc_Q_m(x)) * d.Decimal(0.5)

    def _f_max(self, steps=100):
        result = abs(self._calc_F_m(self.xa))
        for i in range(1, steps):
            val = abs(self._calc_F_m((d.Decimal(i) / d.Decimal(steps)) * self.xa))
            if val > result:
                result = val
        return result

    F_max = property(fget=_f_max)


class Table_13_2(_Table_13_1_and_2_base):
    def __init__(self, xa, v_star, q_3, zero_effective_pressure=False):
        super().__init__(xa, v_star, q_3)
        self.pe = zero_effective_pressure

    def _calc_Z_d_x(self, x):
        num = self._phi_sub_two(self.xa) * self._ber(x) + self._phi_sub_one(self.xa) * self._bei(x)
        den = (self.xa ** 3) * self.Z_a
        return num / den

    def _calc_Z_w_x(self, x):
        num = self._ber_prime(self.xa) * self._ber(x) + self._bei_prime(self.xa) * self._bei(x)
        den = (self.xa ** 2) * self.Z_a
        return num / den

    def _F_t(self, x):
        if self.pe:
            return ((self.xa ** 4) / 2) * ((self._calc_Z_d_x(x)) + self.q_3 * self._calc_Z_w_x(x))
        else:
            return self._calc_Z_w_x(x) * ((self.xa ** 4) / 2)

    def _F_t_min(self, steps=100):
        result = self._F_t(self.xa)
        for i in range(1, steps):
            val = self._F_t((d.Decimal(i) / d.Decimal(steps)) * self.xa)
            if val < result:
                result = val
        return result

    def _F_t_max(self, steps=100):
        result = self._F_t(self.xa)
        for i in range(1, steps):
            val = self._F_t((d.Decimal(i) / d.Decimal(steps)) * self.xa)
            if val > result:
                result = val
        return result

    F_t_min = property(fget=_F_t_min)
    F_t_max = property(fget=_F_t_max)


if __name__ == "__main__":
    calcs = Table_13_1(xa=3.2085,
                       v_star=0.3642,
                       q_3=0.5237018
                       )

    print("***Inputs***")
    print("xa = " + str(calcs.xa))
    print("v_star = " + str(calcs.v_star))
    print("q_3 = " + str(calcs.q_3))
    print("***Outputs***")
    print("Z_a = " + str(calcs.Z_a))
    print("Z_d = " + str(calcs.Z_d))
    print("Z_v = " + str(calcs.Z_v))
    print("Z_w = " + str(calcs.Z_w))
    print("Z_m = " + str(calcs.Z_m))
    print("F_max = " + str(calcs.F_max))
