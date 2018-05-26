
#  TODO: Implement acceptibility tests


class Appendix13_7_bParams:
    def __init__(self,
                 long_side_length_inside,
                 short_side_length_inside,
                 internal_pressure,
                 short_side_thickness,
                 long_side_thickness,
                 long_side_thickness_2,
                 eval_at_outer_walls = False
                 ):
        self.h = short_side_length_inside
        self.H = long_side_length_inside
        self.P = internal_pressure
        self.t_1 = short_side_thickness
        self.t_2 = long_side_thickness
        self.t_22 = long_side_thickness_2
        self.eval_at_outer_walls = eval_at_outer_walls


class Appendix13_7_bCalcs:
    def __init__(self, params: Appendix13_7_bParams):
        self.P = params.P
        self.h = params.h
        self.H = params.H
        self.t_1 = params.t_1
        self.t_2 = params.t_2
        self.t_22 = params.t_22
        self.isOuterWallEval = params.eval_at_outer_walls

    def alpha(self):
        """

        :return: Aspect ratio H / h
        """
        return self.H / self.h

    def I_1(self):
        """

        :return: Moment of inertia of strip thickness t_1
        """
        return (1 / 12.0) * self.t_1 ** 3

    def I_2(self):
        """

        :return: Moment of inertia of strip thickness t_2
        """
        return (1 / 12.0) * self.t_2 ** 3

    def I_22(self):
        """

        :return: Moment of inertia of strip thickness t_22
        """
        return (1 / 12.0) * self.t_22 ** 3

    def k_1(self):
        """

        :return: I_22 / I_2
        """
        return self.I_22() / self.I_2()

    def k_2(self):
        """

        :return: I_22 * alpha / I_1
        """
        return self.I_22() * self.alpha() / self.I_2()

    def K_1(self):
        """

        :return: 2 * k_2 + 3
        """
        return 2.0 * self.k_2() + 3.0

    def K_2(self):
        """

        :return: 3 * k_1 + 2 * k_2
        """
        return 3.0 * self.k_1() + 2.0 * self.k_2()

    def N(self):
        """

        :return: K_1 * K_2 - k_2^2
        """
        return self.K_1() * self.K_2() - self.k_2() ** 2

    def c_1(self):
        """

        :return: The distance from the neutral axis of cross section to extreme fibers. Will return c_i or c_o for its thickness, depending on pressure
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_1

    def c_2(self):
        """

        :return: See c_1
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_2

    def c_22(self):
        """

        :return: see c_1
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_22

    def SmShort(self):
        """

        :return: Short side membrane stress for Figure 13-2(a) Sketch 2 vessels; appendix 13-7 equation 11
        """
        return self.P * self.h / (2.0 * self.t_1)

    def Sm_t_2(self):
        """

        :return: Long side membrane stress for Figure 13-2(a) Sketch 2 vessels of thickness t2; Appendix 13-7 equation 12A
        """
        term_1 = (self.P / (8 * self.N() * self.H * self.t_2))
        constants = (self.K_2() + self.k_2()) - self.k_1() * (self.K_1() + self.k_2()) + (self.alpha() ** 2) * self.k_2() * (self.K_2() - self.K_1())

        return term_1 * (4.0 * self.N() * (self.H ** 2) - 2.0 * (self.h ** 2) * constants)

    def Sm_t_22(self):
        """

        :return: Long side membrane stress for Figure 13-2(a) Sketch 2 vessels of thickness t22; Appendix 13-7 equation 12B
        """
        term_1 = (self.P / (8 * self.N() * self.H * self.t_22))
        constants = -(self.K_2() + self.k_2()) + self.k_1() * (self.K_1() + self.k_2()) - (self.alpha() ** 2) * self.k_2() * (self.K_2() - self.K_1())

        return term_1 * (4.0 * self.N() * (self.H ** 2) - 2.0 * (self.h ** 2) * constants)

    def __Q_bending_geometry_factor(self):
        """

        :return: The factor common to bending stresses at point Q
        """
        return (self.K_2() - self.k_1() * self.k_2()) + (self.alpha() ** 2) * self.k_2() * (self.K_2() - self.k_2())

    def __Q_1_bending_geometry_factor(self):
        """

        :return:  The factor common to bending stresses at point Q_1
        """
        return (self.K_1() * self.k_1() - self.k_2()) + (self.alpha() ** 2) * self.k_2() * (self.K_1() - self.k_2())

    def __short_side_bending_force(self):
        """

        :return: The force factor common to short side bending stresses
        """
        return (self.P * self.c_1() * self.h ** 2) / (4.0 * self.N() * self.I_1())

    def __long_side_bending_force_t_22(self):
        """

        :return: The force factor common to long side bending stresses for thickness t_22
        """
        return (self.P * self.c_22() * self.h ** 2) / (8 * self.N() * self.I_22())

    def __long_side_bending_force_t_2(self):
        """

        :return: The force factor common to long side bending stresses for thickness t_2
        """
        return (self.P * self.c_2() * self.h ** 2) / (8 * self.N() * self.I_2())

    def S_b_Q_short(self):
        """

        :return: Short side bending stress at corner of t_1 and t_22 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 13
        """
        return self.__short_side_bending_force() * self.__Q_bending_geometry_factor()

    def S_b_Q_1_short(self):
        """

        :return: Short side bending stress at corner of t_1 and t_2 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 14
        """
        return self.__short_side_bending_force() * self.__Q_1_bending_geometry_factor()

    def S_b_M(self):
        """

        :return: Long side bending stress at the center of t_22 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 15
        """
        return self.__long_side_bending_force_t_22() * (2 * self.__Q_bending_geometry_factor() - self.N())

    def S_b_M_1(self):
        """

        :return: Long side bending stress at the center of t_2 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 16
        """
        return self.__long_side_bending_force_t_2() * (2 * self.__Q_1_bending_geometry_factor() - self.N())

    def S_b_Q_long(self):
        """

        :return: Long side bending stress at corner of t_1 and t_22 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 17
        """
        return self.__long_side_bending_force_t_22() * self.__Q_bending_geometry_factor()

    def S_b_Q_1_long(self):
        """

        :return: Long side bending stress at corner of t_1 and t_2 for Figure 13-2(a) Sketch 2 vessels; Appendix 13-7 equation 18
        """
        return self.__long_side_bending_force_t_2() * self.__Q_1_bending_geometry_factor()

    def S_T_Q_short(self):
        """

        :return: Total stress at point Q from the short side plate, junction of t_1 and t_22 for Figure 13-2(a) vessels; Appendix 13-7 equation 19
        """
        return self.SmShort() + self.S_b_Q_short()

    def S_T_Q_1_short(self):
        """

        :return: Total stress at point Q1 from the short side plate, junction of t_1 and t_2 for Figure 13-2(a) vessels; Appendix 13-7 equation 20
        """
        return self.SmShort() + self.S_b_Q_1_short()

    def S_T_M(self):
        """

        :return: Total stress at point M, middle of plate of thickness t_22 for Figure 13-2(a) vessels; Appendix 13-7 equation 21
        """
        return self.Sm_t_22() + self.S_b_M()

    def S_T_M_1(self):
        """

        :return: Total stress at point M1, middle of plate of thickness t_2 for Figure 13-2(a) vessels; Appendix 13-7 equation 22
        """
        return self.Sm_t_2() + self.S_b_M_1()

    def S_T_Q_long(self):
        """

        :return: Total stress at point Q from the long side plate, junction of t_1 and t_22 for Figure 13-2(a) vessels; Appendix 13-7 equation 23
        """
        return self.Sm_t_22() + self.S_b_Q_long()

    def S_T_Q_1_long(self):
        """

        :return: Total stress at point Q1 from the long side plate, junction of t_1 and t_2 for Figure 13-2(a) vessels; Appendix 13-7 equation 24
        """
        return self.Sm_t_2() + self.S_b_Q_1_long()

if __name__ == "__main__":
    import copy
    params_inner = Appendix13_7_bParams(
        long_side_length_inside=10,
        short_side_length_inside=5,
        internal_pressure=50,
        short_side_thickness=1,
        long_side_thickness=2,
        long_side_thickness_2=4
    )

    calc_inner = Appendix13_7_bCalcs(params_inner)
    params_outer = copy.deepcopy(params_inner)
    params_outer.eval_at_outer_walls = True
    calc_outer = Appendix13_7_bCalcs(params_outer)

    print("*** Input ***")
    print("P = " + str(calc_inner.P))
    print("h = " + str(calc_inner.h))
    print("H = " + str(calc_inner.H))
    print("t_1 = " + str(calc_inner.t_1))
    print("t_2 = " + str(calc_inner.t_2))
    print("t_22 = " + str(calc_inner.t_22))
    print("\n")
    print("*** Output ***")
    print("")
    print("*** Inner Wall ***")
    print("alpha = " + str(calc_inner.alpha()))
    print("I_1 = " + str(calc_inner.I_1()))
    print("I_2 = " + str(calc_inner.I_2()))
    print("I_22 = " + str(calc_inner.I_22()))
    print("k_1 = " + str(calc_inner.k_1()))
    print("k_2 = " + str(calc_inner.k_2()))
    print("K_1 = " + str(calc_inner.K_1()))
    print("K_2 = " + str(calc_inner.K_2()))
    print("N = " + str(calc_inner.N()))
    print("c_1 = " + str(calc_inner.c_1()))
    print("c_2 = " + str(calc_inner.c_2()))
    print("c_22 = " + str(calc_inner.c_22()))
    print("SmShort = " + str(calc_inner.SmShort()))
    print("Sm_t_2 = " + str(calc_inner.Sm_t_2()))
    print("Sm_t_22 = " + str(calc_inner.Sm_t_22()))
    print("S_b_Q_short = " + str(calc_inner.S_b_Q_short()))
    print("S_b_Q_1_short = " + str(calc_inner.S_b_Q_1_short()))
    print("S_b_M = " + str(calc_inner.S_b_M()))
    print("S_b_M_1 = " + str(calc_inner.S_b_M_1()))
    print("S_b_Q_long = " + str(calc_inner.S_b_Q_long()))
    print("S_b_Q_1_long = " + str(calc_inner.S_b_Q_1_long()))
    print("S_T_Q_short = " + str(calc_inner.S_T_Q_short()))
    print("S_T_Q_1_short = " + str(calc_inner.S_T_Q_1_short()))
    print("S_T_M = " + str(calc_inner.S_T_M()))
    print("S_T_M_1 = " + str(calc_inner.S_T_M_1()))
    print("S_T_Q_long = " + str(calc_inner.S_T_Q_long()))
    print("S_T_Q_1_long = " + str(calc_inner.S_T_Q_1_long()))
    print("")
    print("*** Outer Wall ***")
    print("alpha = " + str(calc_outer.alpha()))
    print("I_1 = " + str(calc_outer.I_1()))
    print("I_2 = " + str(calc_outer.I_2()))
    print("I_22 = " + str(calc_outer.I_22()))
    print("k_1 = " + str(calc_outer.k_1()))
    print("k_2 = " + str(calc_outer.k_2()))
    print("K_1 = " + str(calc_outer.K_1()))
    print("K_2 = " + str(calc_outer.K_2()))
    print("N = " + str(calc_outer.N()))
    print("c_1 = " + str(calc_outer.c_1()))
    print("c_2 = " + str(calc_outer.c_2()))
    print("c_22 = " + str(calc_outer.c_22()))
    print("SmShort = " + str(calc_outer.SmShort()))
    print("Sm_t_2 = " + str(calc_outer.Sm_t_2()))
    print("Sm_t_22 = " + str(calc_outer.Sm_t_22()))
    print("S_b_Q_short = " + str(calc_outer.S_b_Q_short()))
    print("S_b_Q_1_short = " + str(calc_outer.S_b_Q_1_short()))
    print("S_b_M = " + str(calc_outer.S_b_M()))
    print("S_b_M_1 = " + str(calc_outer.S_b_M_1()))
    print("S_b_Q_long = " + str(calc_outer.S_b_Q_long()))
    print("S_b_Q_1_long = " + str(calc_outer.S_b_Q_1_long()))
    print("S_T_Q_short = " + str(calc_outer.S_T_Q_short()))
    print("S_T_Q_1_short = " + str(calc_outer.S_T_Q_1_short()))
    print("S_T_M = " + str(calc_outer.S_T_M()))
    print("S_T_M_1 = " + str(calc_outer.S_T_M_1()))
    print("S_T_Q_long = " + str(calc_outer.S_T_Q_long()))
    print("S_T_Q_1_long = " + str(calc_outer.S_T_Q_1_long()))