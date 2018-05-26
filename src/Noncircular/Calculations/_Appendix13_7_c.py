import math

#  TODO: Implement acceptibility tests

class Appendix13_7_cParams:
    def __init__(
            self,
            internal_pressure,
            corner_radius,
            short_side_half_length,
            long_side_half_length,
            thickness,
            eval_at_outer_walls = False):

        self.P = internal_pressure
        self.R = corner_radius
        self.L_1 = short_side_half_length
        self.L_2 = long_side_half_length
        self.t_1 = thickness
        self.eval_at_outer_walls = eval_at_outer_walls


class Appendix13_7_cCalcs:
    def __init__(self, params: Appendix13_7_cParams):
        self.P = params.P
        self.R = params.R
        self.L_1 = params.L_1
        self.L_2 = params.L_2
        self.t_1 = params.t_1
        self.isOuterWallEval = params.eval_at_outer_walls

    def c(self):
        """

        :return: The distance from the neutral axis of cross section to extreme fibers. Will return c_i or c_o for its thickness, depending on pressure
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_1

    def I_1(self):
        return (1 / 12.0) * self.t_1 ** 3

    def alpha3(self):
        return self.L_2 / self.L_1

    def phi(self):
        return self.R / self.L_1

    def K_3(self):
        """

        :return: Equation 40
        """
        return (-1.0) * (self.L_1 ** 2) * (
            6.0 * (self.phi() ** 2) * self.alpha3()
            - 3.0 * math.pi * (self.phi() ** 2)
            + 6.0 * (self.phi() ** 2)
            + (self.alpha3() ** 3)
            + (3.0 * self.alpha3() ** 2)
            - 6.0 * self.phi()
            - 2.0
            + 1.5 * math.pi * self.phi() * (self.alpha3() ** 2)
            + 6.0 * self.phi() * self.alpha3()
        ) / (3.0 * (2.0 * self.alpha3() + math.pi * self.phi() + 2.0))

    def M_A(self):
        """

        :return: Equation 38
        """
        return self.P * self.K_3()

    def M_r(self):
        """

        :return: equation 39
        """
        raise ValueError("Looks like it's time to implement M_r")

    def S_m_C(self):
        """

        :return: Short side membrane stress at point C for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 25
        """
        return (self.P * (self.R + self.L_2)) / self.t_1

    def S_m_D(self):
        """

        :return: Same as S_m_C
        """
        return self.S_m_C()

    def S_m_A(self):
        """

        :return: Long side membrane stress at point A for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 26
        """
        return (self.P *(self.L_1 + self.R)) / self.t_1

    def S_m_B(self):
        """

        :return: Same as S_m_A
        """
        return self.S_m_A()

    def S_m_BC(self):
        """

        :return: Membrane stress in radius, between points B and C for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 27
        """
        return (self.P / self.t_1) * (math.sqrt((self.L_2 ** 2) + self.L_1 ** 2) + self.R)

    def S_b_C(self):
        """

        :return: Bending stress at C for short side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 28
        """
        return (self.c() / (2.0 * self.I_1())) * (2.0 * self.M_A() + self.P * (2 * self.R * self.L_2 - 2.0 * self.R * self.L_1 + self.L_2 ** 2))

    def S_b_D(self):
        """

        :return: Bending stress at D for short side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 29
        """
        return (self.c() / (2.0 * self.I_1())) * (2.0 * self.M_A() + self.P * ((self.L_2 ** 2) + 2 * self.R * self.L_2 - 2.0 * self.R * self.L_1 + self.L_2 ** 2))

    def S_b_A(self):
        """

        :return: Bending stress at point A for long side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 30
        """
        return self.M_A() * self.c() / self.I_1()

    def S_b_B(self):
        """

        :return: Bending stress at point B for long side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 31
        """
        return (self.c() / (2 * self.I_1())) * (2 * self.M_A() + self.P * self.L_2 ** 2)

    def S_b_BC(self):
        """

        :return: Max bending stress between points B and C for corner sections for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 32
        """
        maxStressTheta = math.atan(self.L_1 / self.L_2)
        geom = self.c() / self.I_1()
        moment = 0.5 * (2 * self.M_A() + self.P * (2 * self.R * (self.L_2 * math.cos(maxStressTheta) - self.L_1 * (1 - math.sin(maxStressTheta))) + self.L_2 ** 2))

        return geom * moment

    def S_T_C(self):
        """

        :return: Total stress at point C for short side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 33
        """
        return self.S_m_C() + self.S_b_C()

    def S_T_D(self):
        """

        :return: Total stress at point D for short side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 34
        """
        return self.S_m_D() + self.S_b_D()

    def S_T_A(self):
        """

        :return: Total stress at point A for long side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 35
        """
        return self.S_m_A() + self.S_b_A()

    def S_T_B(self):
        """

        :return: Total stress at point B for long side plate for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 36
        """
        return self.S_m_B() + self.S_b_B()

    def S_T_BC(self):
        """

        :return: Total stress between points B and C for corner sections for Figure 13-2(a) Sketch 3 vessels; appendix 13-7 equation 37
        """
        return self.S_m_BC() + self.S_b_BC()

if __name__ == "__main__":
    import copy
    params_inner = Appendix13_7_cParams(
        internal_pressure=100,
        corner_radius=3,
        short_side_half_length=5,
        long_side_half_length=10,
        thickness=1
    )

    calc_inner = Appendix13_7_cCalcs(params_inner)

    params_outer = copy.deepcopy(params_inner)
    params_outer.eval_at_outer_walls = True
    calc_outer = Appendix13_7_cCalcs(params_outer)


    print("*** Input ***")
    print("P = " + str(params_inner.P))
    print("R = " + str(params_inner.R))
    print("L_1 = " + str(params_inner.L_1))
    print("L_2 = " + str(params_inner.L_2))
    print("t_1 = " + str(params_inner.t_1))
    print("")

    print("*** Output ***")
    print("")
    print("*** Inner Walls ***")
    print("c = " + str(calc_inner.c()))
    print("I_1 = " + str(calc_inner.I_1()))
    print("alpha3 = " + str(calc_inner.alpha3()))
    print("phi = " + str(calc_inner.phi()))
    print("K_3 = " + str(calc_inner.K_3()))
    print("M_A = " + str(calc_inner.M_A()))
    # print("M_r = " + str(calc_inner.M_r()))
    print("S_m_C = " + str(calc_inner.S_m_C()))
    print("S_m_D = " + str(calc_inner.S_m_D()))
    print("S_m_A = " + str(calc_inner.S_m_A()))
    print("S_m_B = " + str(calc_inner.S_m_B()))
    print("S_m_BC = " + str(calc_inner.S_m_BC()))
    print("S_b_C = " + str(calc_inner.S_b_C()))
    print("S_b_D = " + str(calc_inner.S_b_D()))
    print("S_b_A = " + str(calc_inner.S_b_A()))
    print("S_b_B = " + str(calc_inner.S_b_B()))
    print("S_b_BC = " + str(calc_inner.S_b_BC()))
    print("S_T_C = " + str(calc_inner.S_T_C()))
    print("S_T_D = " + str(calc_inner.S_T_D()))
    print("S_T_A = " + str(calc_inner.S_T_A()))
    print("S_T_B = " + str(calc_inner.S_T_B()))
    print("S_T_BC = " + str(calc_inner.S_T_BC()))
    print("")
    print("*** Outer Walls ***")
    print("c = " + str(calc_outer.c()))
    print("I_1 = " + str(calc_outer.I_1()))
    print("alpha3 = " + str(calc_outer.alpha3()))
    print("phi = " + str(calc_outer.phi()))
    print("K_3 = " + str(calc_outer.K_3()))
    print("M_A = " + str(calc_outer.M_A()))
    # print("M_r = " + str(calc_outer.M_r()))
    print("S_m_C = " + str(calc_outer.S_m_C()))
    print("S_m_D = " + str(calc_outer.S_m_D()))
    print("S_m_A = " + str(calc_outer.S_m_A()))
    print("S_m_B = " + str(calc_outer.S_m_B()))
    print("S_m_BC = " + str(calc_outer.S_m_BC()))
    print("S_b_C = " + str(calc_outer.S_b_C()))
    print("S_b_D = " + str(calc_outer.S_b_D()))
    print("S_b_A = " + str(calc_outer.S_b_A()))
    print("S_b_B = " + str(calc_outer.S_b_B()))
    print("S_b_BC = " + str(calc_outer.S_b_BC()))
    print("S_T_C = " + str(calc_outer.S_T_C()))
    print("S_T_D = " + str(calc_outer.S_T_D()))
    print("S_T_A = " + str(calc_outer.S_T_A()))
    print("S_T_B = " + str(calc_outer.S_T_B()))
    print("S_T_BC = " + str(calc_outer.S_T_BC()))