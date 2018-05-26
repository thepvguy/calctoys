

class Appendix13_9_bParams:
    def __init__(self,
                 dist_from_short_side_to_stay_plate,
                 short_side_length,
                 internal_pressure,
                 short_side_thickness,
                 long_side_thickness,
                 stay_plate_thickness,
                 eval_at_outer_walls = False
                 ):
        self.h = dist_from_short_side_to_stay_plate
        self.H = short_side_length
        self.P = internal_pressure
        self.t_1 = short_side_thickness
        self.t_2 = long_side_thickness
        self.t_3 = stay_plate_thickness
        self.eval_at_outer_walls = eval_at_outer_walls

class Appendix13_9_bCalcs:
    def __init__(self, params):
        self.h = params.h
        self.H = params.H
        self.P = params.P
        self.t_1 = params.t_1
        self.t_2 = params.t_2
        self.t_3 = params.t_3
        self.isOuterWallEval = params.eval_at_outer_walls

    def I_1(self):
        """

        :return: longitudinal unit moment of inertia of side with thickness t1
        """
        return (1 / 12.0) * self.t_1 ** 3

    def I_2(self):
        """

        :return: longitudinal unit moment of inertia of side with thickness t2
        """
        return (1 / 12.0) * self.t_2 ** 3

    def alpha(self):
        """

        :return: Rectangular vessel parameter
        """
        return self.H / self.h

    def K(self):
        """

        :return: Vessel parameter
        """
        return (self.I_2() / self.I_1()) * self.alpha()

    def c_1(self):
        """

        :return: The distance from the neutral axis to the extreme fibers of thickness t1
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_1

    def c_2(self):
        """

        :return: The distance from the neutral axis to the extreme fibers of thickness t2
        """
        sign = 1
        if self.isOuterWallEval:
            sign = -1
        return 0.5 * sign * self.t_2

    def __S_m_common(self):
        """

        :return: Common factor between short side and stay plate membrane stresses
        """

        return self.P * self.h * ((2 + self.K() * (5 - self.alpha() ** 2)) / (1 + 2 * self.K()))

    def SmShort(self):
        """

        :return: Short Side Membrane stress for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 1
        """
        return self.__S_m_common() / (4 * self.t_1)

    def SmLong(self):
        """

        :return: Long Side Membrane stress for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 2
        """
        return self.P * self.H / (2 * self.t_2)

    def SmStay(self):
        """

        :return: Stay Plate Membrane stress for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 3
        """
        return self.__S_m_common() / (2 * self.t_3)

    def __S_b_common(self):
        """

        :return: Factor common to bending stresses
        """
        return (self.h ** 2) * ((1 + 2 * (self.alpha() ** 2) * self.K()) / (1 + 2 * self.K()))

    def S_b_N(self):
        """

        :return: Short side plate bending stress at midpoint for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 4
        """
        return ((self.P * self.c_1()) / (24 * self. I_1())) * ((-3 * self.H ** 2) + 2 * self.__S_b_common())

    def S_b_Q_short(self):
        """

        :return: Short side plate bending stress at corner for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 5
        """
        return ((self.P * self.c_1()) / (12 * self.I_1())) * self.__S_b_common()

    def S_b_M(self):
        """

        :return: Long side plate bending stress at midspan for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 6
        """
        return ((self.P * self.c_2() * self.h ** 2) / (12 * self.I_2())) * ((1 + self.K() * (3 - self.alpha() ** 2)) / (1 + 2 * self.K()))

    def S_b_Q_long(self):
        """

        :return: Long side plate bending stress at corner for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 7
        """
        return ((self.P * self.c_2()) / (12 * self.I_2())) * self.__S_b_common()

    def S_T_N(self):
        """

        :return: Short side plate total stress at midpoint for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 8
        """
        return self.SmShort() + self.S_b_M()

    def S_T_Q_short(self):
        """

        :return: Short side plate total stress at corner for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 9
        """
        return self.SmShort() + self.S_b_Q_short()

    def S_T_M(self):
        """

        :return: Long side plate total stress at midpoint for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 10
        """
        return self.SmLong() + self.S_b_M()

    def S_T_Q_long(self):
        """

        :return: Short side plate total stress at corner for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 11
        """
        return self.SmLong() + self.S_b_Q_long()

    def S_T_stay(self):
        """

        :return: Stay plate total stress for Figure 13-2(a) Sketch 7 vessels; Appendix 13-9 equation 12
        """
        return self.SmStay()

if __name__ == "__main__":
    import copy
    params_inner = Appendix13_9_bParams(
        dist_from_short_side_to_stay_plate=5,
        short_side_length=5,
        internal_pressure=100,
        short_side_thickness=1,
        long_side_thickness=2,
        stay_plate_thickness=0.5
    )
    params_outer = copy.deepcopy(params_inner)
    params_outer.evalAtOuterWalls = True

    calc_inner = Appendix13_9_bCalcs(params_inner)
    calc_outer = Appendix13_9_bCalcs(params_outer)

    print("*** Evaluated at inner walls ***")
    print("c_1 = " + str(calc_inner.c_1()))
    print("c_2 = " + str(calc_inner.c_2()))
    print("alpha = " + str(calc_inner.alpha()))
    print("I_1 = " + str(calc_inner.I_1()))
    print("I_2 = " + str(calc_inner.I_2()))
    print("K = " + str(calc_inner.K()))
    print("SmShort = " + str(calc_inner.SmShort()))
    print("SmLong = " + str(calc_inner.SmLong()))
    print("S_b_N = " + str(calc_inner.S_b_N()))
    print("S_b_Q_short = " + str(calc_inner.S_b_Q_short()))
    print("S_b_M = " + str(calc_inner.S_b_M()))
    print("S_b_Q_long = " + str(calc_inner.S_b_Q_long()))
    print("S_T_N = " + str(calc_inner.S_T_N()))
    print("S_T_Q_short = " + str(calc_inner.S_T_Q_short()))
    print("S_T_M = " + str(calc_inner.S_T_M()))
    print("S_T_Q_long = " + str(calc_inner.S_T_Q_long()))
    print("S_T_stay = " + str(calc_inner.S_T_stay()))
    print("")
    print("*** Evaluated at outer walls ***")
    print("c_1 = " + str(calc_outer.c_1()))
    print("c_2 = " + str(calc_outer.c_2()))
    print("alpha = " + str(calc_outer.alpha()))
    print("I_1 = " + str(calc_outer.I_1()))
    print("I_2 = " + str(calc_outer.I_2()))
    print("K = " + str(calc_outer.K()))
    print("SmShort = " + str(calc_outer.SmShort()))
    print("SmLong = " + str(calc_outer.SmLong()))
    print("S_b_N = " + str(calc_outer.S_b_N()))
    print("S_b_Q_short = " + str(calc_outer.S_b_Q_short()))
    print("S_b_M = " + str(calc_outer.S_b_M()))
    print("S_b_Q_long = " + str(calc_outer.S_b_Q_long()))
    print("S_T_N = " + str(calc_outer.S_T_N()))
    print("S_T_Q_short = " + str(calc_outer.S_T_Q_short()))
    print("S_T_M = " + str(calc_outer.S_T_M()))
    print("S_T_Q_long = " + str(calc_outer.S_T_Q_long()))
    print("S_T_stay = " + str(calc_outer.S_T_stay()))