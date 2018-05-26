"""
    Un-reinforced vessels of rectangular cross section
"""

# from . import _Appendix13Common
from ._Appendix13_6 import MultiDiameterHole
#  TODO: Implement acceptibility tests
#  TODO: Include ligament efficiencies

class Appendix13_7_aParams:
    def __init__(self,
                 long_side_length_inside,
                 short_side_length_inside,
                 internal_pressure,
                 short_side_thickness,
                 long_side_thickness,
                 allowable_stress,
                 joint_efficiency,
                 eval_at_outer_walls=False
                 ):

        self.h = long_side_length_inside
        self.H = short_side_length_inside
        self.P = internal_pressure
        self.t_1 = short_side_thickness
        self.t_2 = long_side_thickness
        self.S = allowable_stress
        self.E = joint_efficiency
        self.evalAtOuterWalls = eval_at_outer_walls


class Appendix13_7_aCalcs:
    def __init__(self, params: Appendix13_7_aParams):
        self.h = params.h
        self.H = params.H
        self.P = params.P
        self.t_1 = params.t_1
        self.t_2 = params.t_2
        self.S = params.S
        self.E = params.E
        self.isOuterWallEval = params.evalAtOuterWalls

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

    def SmShort(self):
        """

        :return: Short side membrane stress for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 1
        """
        return self.P * self.h / (2.0 * self.t_1)

    def SmLong(self):
        """

        :return: Long side membrane stress for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 2
        """
        return self.P * self.H / (2.0 * self.t_2)

    def S_b_N(self):
        """

        :return: Short side bending stress for Figure 13-2(a) Sketch 1 vessels in midspan; appendix 13-7 equation 3
        """
        return (self.P * self.c_1() / (12 * self.I_1())) * (
        -1.5 * (self.H ** 2) + (self.h ** 2) * ((1 + (self.alpha() ** 2) * self.K()) / (1 + self.K())))

    def S_b_Q_short(self):
        """

        :return: Short side bending stress for Figure 13-2(a) Sketch 1 vessels at corner; appendix 13-7 equation 4
        """
        return ((self.P * (self.h ** 2) * self.c_1()) / (12 * self.I_1())) * (
        (1 + (self.alpha() ** 2) * self.K()) / (1 + self.K()))

    def S_b_M(self):
        """

        :return: Long side bending stress for Figure 13-2(a) Sketch 1 vessels in midspan; appendix 13-7 equation 5
        """
        return ((self.P * (self.h ** 2) * self.c_2()) / (12 * self.I_2())) * (
        -1.5 + ((1 + (self.alpha() ** 2) * self.K()) / (1 + self.K())))

    def S_b_Q_long(self):
        """

        :return: Long side bending stress for Figure 13-2(a) Sketch 1 vessels at corner; appendix 13-7 equation 6
        """
        return ((self.P * (self.h ** 2) * self.c_2()) / (12 * self.I_2())) * (
        (1 + (self.alpha() ** 2) * self.K()) / (1 + self.K()))

    def S_T_N(self):
        """

        :return: Total stress at midspan of short side of vessel for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 7
        """
        return self.SmShort() + self.S_b_N()

    def S_T_Q_short(self):
        """

        :return: Total stress at corner of short side of vessel for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 8
        """
        return self.SmShort() + self.S_b_Q_short()

    def S_T_M(self):
        """

        :return: Total stress at midspan of long side of vessel for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 9
        """
        return self.SmLong() + self.S_b_M()

    def S_T_Q_long(self):
        """

        :return: Total stress at corner of long side of vessel for Figure 13-2(a) Sketch 1 vessels; appendix 13-7 equation 10
        """
        return self.SmLong() + self.S_b_Q_long()

    def endStress(self):
        return self.h * self.H * self.P



if __name__ == "__main__":
    import copy
    params_inner = Appendix13_7_aParams(
        long_side_length_inside=9.5,
        short_side_length_inside=7.375,
        internal_pressure=400,
        short_side_thickness=0.875,
        long_side_thickness=0.875,
        allowable_stress=20000,
        joint_efficiency=1
    )
    params_outer = copy.deepcopy(params_inner)
    params_outer.evalAtOuterWalls = True

    calc_inner = Appendix13_7_aCalcs(params_inner)
    calc_outer = Appendix13_7_aCalcs(params_outer)

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
