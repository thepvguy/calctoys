import math
import enum
from typing import Optional


class Appendix2Params:
    def __init__(
            self,
            other
    ):
        pass


class Appendix2FlangeCalcs:
    """
        This is intended to be the flange base class.
        ***Incomplete.***
        Refer to ASME BPVC VIII-1 Appendix 2 for a list of variables -
        anything that is calculated is a property of this object

        The functions in this class roughly mirror the variable names in Appendix 2.
        In cases where there is both an operating and bolting case, the word "operating" and "bolting" is appended
        to the end of the function name, respectively.
        """

    # TODO: External pressure
    # TODO: Reverse flanges
    # TODO: Split flanges
    # TODO: Allowable Stresses
    # TODO: Split flanges
    # TODO: Noncircular bolt pattern
    # TODO: Nut Stops
    # TODO: Material checks

    def __init__(
            self,
            params):

        # Design conditions
        self.P = params.internal_pressure
        self.S_a = params.bolt_ambient_allowable_stress
        self.S_b = params.bolt_operating_allowable_stress
        self.S_f = params.flange_operating_allowable_stress
        self.S_n = params.nozzle_neck_operating_allowable
        self.E_a = params.flange_ambient_modulus_of_elasticity
        self.E_f = params.flange_operating_modulus_of_elasticity
        self._custom_W_m_1_active = params.custom_Wm1_active
        self.custom_W_m_1 = params.custom_Wm1
        self._custom_W_m_2_active = params.custom_Wm2_active
        self.custom_W_m_2 = params.custom_Wm2

        # Gasket stuff
        self.G_od = params.gasket_OD
        self.G_id = params.gasket_ID
        self.gasket_m = params.gasket_m
        self.gasket_y = params.gasket_y
        self.T_gasket = params.gasket_thickness
        self.gasket_row = params.facing_sketch
        self.gasket_col = params.facing_column
        self.w = params.w
        self.nubbin_height = params.nubbin_height
        self.raised_face_type = params.raised_face_type
        self.facing_dia = params.rf_dia

        # Flange dimensions
        self.A = params.flange_OD
        self.B = params.flange_ID
        self.t = params.flange_thickness
        self.numBolts = params.num_bolts
        self.C = params.bolt_circle_diameter
        self.hasHub = params.has_hub
        if self.hasHub:
            self.hub = params.hub
        else:
            self.hub = None
        self.a = params.bolt_dismeter
        self.boltRootArea = params.bolt_root_area

        self.attachment_sketch = params.attachment_sketch
        self.__custom_rigidity_factor = params.custom_rigidity
        self.units = params.units

    @property
    def C_b(self):
        if self.units == Units.Imperial:
            return 0.5
        elif self.units == Units.Metric:
            return 2.5
        else:
            raise ValueError("Improper value %r for property 'units'" % self.units)

    @property
    def b(self):
        """
        The effective gasket or joint-contact-surface seating width
        :return:
        """
        if self.b_o > .25:
            return self.C_b * math.sqrt(self.b_o)
        elif self.b_o <= .25:
            return self.b_o
        else:
            raise ValueError("Invalid value of b nought")

    @property
    def N(self):
        """

        :return: Width used to determine the basic gasket seating width b_o,
        based upon the possible contact width of the gasket
        """
        #  Not entirely correct, but should work OK for now
        return (self.G_od - self.G_id) / 2  # TODO : account for gasket falling off of facing

    @property
    def B_1(self):
        if self.attachment_sketch in loose_flanges or self.f_actual < 1:
            return self.B + self.hub.g_1
        elif (self.attachment_sketch in integral_flanges or optional_flanges) and (self.f_actual >= 1):
            return self.B + self.hub.g_1
        else:
            raise ValueError("Invalid value %r for property 'attachment sketch'" % self.attachment_sketch)

    @property
    def b_o(self):
        """
        Basic gasket seating width
        :return:
        """
        #  from table 2-5.2
        if self.gasket_col == 1:
            if self.gasket_row == Table_2_5_2_Sketch.s_1a or self.gasket_row == Table_2_5_2_Sketch.s_1b:
                return self.N / 2
            elif self.gasket_row == Table_2_5_2_Sketch.s_1c or self.gasket_row == Table_2_5_2_Sketch.s_1d:
                retval = (self.w + self.T_gasket) / 2
                maximum = (self.w + self.N) / 4

                if retval >= maximum:
                    return maximum
                elif retval < maximum:
                    return retval
                else:
                    raise ValueError("%r, %r, or %r is not a valid value for properties w, N, or T" %
                                     (self.w, self.N, self.T_gasket))
            elif self.gasket_row == Table_2_5_2_Sketch.s_2:
                return (self.w + self.N) / 4.0

            elif self.gasket_row == Table_2_5_2_Sketch.s_3:
                return self.N / 4.0

            elif self.gasket_row == Table_2_5_2_Sketch.s_4:
                return self.N * (3.0 / 8.0)

            elif self.gasket_row == Table_2_5_2_Sketch.s_5:
                return self.N / 4.0

            elif self.gasket_row == Table_2_5_2_Sketch.s_6:
                return self.w / 8.0
            else:
                raise ValueError("%r is not an appropriate value for property 'Sketch Row'" % self.gasket_row)

        elif self.gasket_col == 2:
            if self.gasket_row == Table_2_5_2_Sketch.s_1a or self.gasket_row == Table_2_5_2_Sketch.s_1b:
                return self.N / 2
            elif self.gasket_row == Table_2_5_2_Sketch.s_1c or self.gasket_row == Table_2_5_2_Sketch.s_1d:
                retval = (self.w + self.T_gasket) / 2
                _max = (self.w + self.N) / 4

                if retval >= _max:
                    return _max
                elif retval < _max:
                    return retval
                else:
                    raise ValueError("%r, %r, or %r is not a valid value for properties w, N, or T" %
                                     (self.w, self.N, self.T_gasket))

            elif self.gasket_row == Table_2_5_2_Sketch.s_2:
                return (self.w + 3.0 * self.N) / 8.0

            elif self.gasket_row == Table_2_5_2_Sketch.s_3:
                return 3.0 * self.N / 8.0

            elif self.gasket_row == Table_2_5_2_Sketch.s_4:
                return self.N * (7.0 / 16.0)

            elif self.gasket_row == Table_2_5_2_Sketch.s_5:
                return 3.0 * self.N / 8.0

            else:
                raise ValueError("%r is not an appropriate value for property 'Sketch Row'" % self.gasket_row)

        else:
            raise ValueError("%r is not an appropriate value of property 'Sketch Column'" % self.gasket_col)

    def effective_gasket_od(self):
        if self.raised_face_type == FacingType.no_facing:
            return self.G_od
        if self.raised_face_type in [FacingType.confined_face, FacingType.raised_face] :
            return min(self.facing_dia, self.G_od)
            # TODO: warn here

    def effective_gasket_id(self):
        return max(self.G_id, self.B)

    @property
    def G(self):
        """
        Diameter at location of gasket load reaction.
        :return:
        """
        if self.b_o > .25:
            return self.G_od - 2 * self.b
        elif self.b_o <= .25:
            return (self.G_od - self.G_id) / 2.0
        else:
            raise ValueError("Invalid value %r for property 'b_o'" % self.b_o)

    @property
    def H(self):
        """
        Total hydrostatic end force
        :return:
        """
        return 0.785 * (self.G ** 2) * self.P

    @property
    def H_p(self):
        """
        Total joint-contact surface compression load
        :return:
        """
        return 2 * self.b * 3.14 * self.G * self.gasket_m * self.P

    @property
    def W_m_1(self):
        """
        Minimum required bolt load for the operating conditions.
        :return:
        """
        if self._custom_W_m_1_active is False:

            return self.H + self.H_p
        elif self._custom_W_m_1_active is True:
            return self.custom_W_m_1
        else:
            raise ValueError("Invalid flag for _custom_W_m_1")

    @W_m_1.setter
    def W_m_1(self, value):
        """
        Sets a custom value for W_m_1. Code use is for tubesheets gasketed both sides.
        :param value: Double
        :return:
        """
        self._custom_W_m_1_active = True
        self.custom_W_m_1 = value

    @property
    def W_m_2(self):
        """
        Minimum requird bolt load for gasket seating.
        :return:
        """
        if self._custom_W_m_2_active is False:
            return 3.14 * self.b * self.G * self.gasket_y
        elif self._custom_W_m_2_active is True:
            return self.custom_W_m_2
        else:
            raise ValueError("Invalid flag for parameter _custom_W_m_2_active")

    @W_m_2.setter
    def W_m_2(self, value):
        """
        Sets a custom value for W_m_2. Code use is for tubesheets gasketed both sides.
        :param value:
        :return:
        """

        self._custom_W_m_2_active = True
        self.custom_W_m_2 = value

    @property
    def B_s_max(self):
        """
        Max bolt spacing
        :return:
        """
        return 2 * self.a + ((6 * self.t) / (self.gasket_m + 0.5))

    @property
    def W_operating(self):
        """
        Formula from 2-4(e)(4) flange design bolt load, operating
        :return:
        """
        return self.W_m_1

    @property
    def W_seating(self):
        """
        Formula from 2-4(e)(5), flange design bolt load, seating
        :return:
        """
        return (self.A_m + self.A_b) * self.S_a * 0.5

    @property
    def A_m(self):
        """
        Total required bolt area
        :return:
        """
        return max(self.A_m_1, self.A_m_2)

    @property
    def A_b(self):
        """
        Bolt root area
        :return:
        """
        return self.boltRootArea * self.numBolts

    @property
    def A_m_1(self):
        """
        Total cross sectional area of bolts at root of thread or section with least diameter required for operating
        :return:
        """
        return self.W_m_1 / self.S_a

    @property
    def A_m_2(self):
        """
        Total cross sectional area of bolts at root of thread or section with least diameter required for gasket seating
        :return:
        """
        return self.W_m_2 / self.S_b

    @property
    def M_o_seating(self):
        """
        Total moment in the seating condition
        :return:
        """
        return self.W_seating * 0.5 * (self.C * self.G)

    @property
    def M_o_operating(self):
        """
        Total moment in the operating condition
        :return:
        """
        return self.M_D + self.M_G_operating + self.M_T

    @property
    def H_D(self):
        """
        Hydrostatic end force on the inside of the flange
        :return:
        """
        return 0.785 * (self.B ** 2) * self.P

    @property
    def W_(self):
        return max(self.W_operating, self.W_seating)

    @property
    def H_G_operating(self):
        """
        Gasket load, operating
        :return:
        """
        return self.W_operating - self.H

    @property
    def H_G_seating(self):
        """
        Gasket load, seating
        :return:
        """
        return self.W_seating - self.H

    @property
    def H_T(self):
        """
        Difference between total hydrostatic end force and hydrostatic end force on inside of flange
        :return:
        """
        return self.H - self.H_D

    @property
    def R(self):
        return 0.5 * (self.C - self.B) - self.hub.g_1

    @property
    def h_D(self):
        """
        Radial distance between bolt circle and H_d radius
        :return:
        """
        if self.attachment_sketch in integral_flanges:
            return self.R + 0.5 * self.hub.g_1
        elif self.attachment_sketch in loose_flanges or \
                        self.attachment_sketch in optional_flanges:
            return 0.5 * (self.C - self.B)
        else:
            raise ValueError("Invalid value '%r' for attribute 'attachment_sketch'" % self.attachment_sketch)

    @property
    def h_G(self):
        """
        Radial distance between gasket load reaction and bolt circle
        :return:
        """
        return (self.C - self.G) * 0.5

    @property
    def h_T(self):
        """
        radial distance between bot circle and H_t
        :return:
        """
        if self.attachment_sketch in integral_flanges:
            return 0.5 * (self.R + self.hub.g_1 + self.h_G)
        elif self.attachment_sketch in loose_flanges or \
                        self.attachment_sketch in optional_flanges:

            if self.attachment_sketch in [Table_2_5_2_Sketch.sketch_1, Table_2_5_2_Sketch.sketch_1a]:
                return 0.5 * (self.C - self.G)
            else:
                return 0.5 * (self.h_D + self.h_G)
        else:
            raise ValueError("Invalid value '%r' for attribute 'attachment_sketch'" % self.attachment_sketch)

    @property
    def M_D(self):
        """
        Moment due to hydrostatic end force
        :return:
        """
        return self.H_D * self.h_D

    @property
    def M_T(self):
        """
        Moment due to difference between hydrostatic end force and HEF on inside of flange
        :return:
        """
        return self.H_T * self.h_T

    @property
    def M_G_seating(self):
        """
        Moment due to gasket load, seating
        :return:
        """
        return self.H_G_seating * self.h_G

    @property
    def M_G_operating(self):
        """
        Moment due to gasket load, operating
        :return:
        """
        return self.H_G_operating * self.h_G

    @property
    def B_S(self):
        return self.C * math.pi / self.numBolts

    @property
    def B_SC(self):
        """
        Bolt spacing correction factor
        :return:
        """
        return math.sqrt(self.B_S / (2 * self.a + self.t))

    @property
    def V(self):
        assert (self.hub is not None)
        if self.hub.g_1 == self.hub.g_o:
            return 0.550103
        else:
            table = Table_2_7_1(self)
            return table.E_4 / (((2.73 / table.C) ** 0.25) * ((1 + table.A) ** 3))

    @property
    def V_L(self):
        assert (self.hub is not None)
        table = Table_2_7_1(self)
        return (0.25 - 0.2 * table.C_24 - 1.5 * table.C_21 - table.C_18) / (
            ((2.73 / table.C) ** 0.25) * (1 + table.A) ** 3)

    @property
    def h_o(self):
        return math.sqrt(self.B * self.hub.g_o)

    @property
    def d(self):
        """
        A factor
        :return:
        """
        if self.attachment_sketch in integral_flanges or \
                        self.attachment_sketch in optional_flanges:
            return (self.U / self.V) * self.h_o * (self.hub.g_o ** 2)
        elif self.attachment_sketch in loose_flanges:
            return (self.U / self.V_L) * self.h_o * (self.hub.g_o ** 2)
        else:
            raise ValueError("Inappropriate value %r for attribute 'Attachment Sketch'" % self.attachment_sketch)

    @property
    def F(self):
        assert (self.hub is not None)
        if self.hub.g_1 == self.hub.g_o:
            return 0.908920
        else:
            table = Table_2_7_1(self)
            return table.E_6 / (((table.C / 2.73) ** 0.25) * (((1 + table.A) ** 3) / table.C))

    @property
    def F_L(self):
        assert (self.hub is not None)
        table = Table_2_7_1(self)
        num1 = table.C_18 * (0.5 + (table.A / 6))
        num2 = table.C_21 * (0.25 + (11 * table.A / 84))
        num3 = table.C_24 * ((1 / 70) + (table.A / 105))
        num4 = (1 / 40) + (table.A / 72)
        den = ((table.C / 2.73) ** 0.25) * (((1 + table.A) ** 3) / table.C)
        return (num1 + num2 + num3 - num4) / den

    @property
    def e(self):
        """
        A factor
        :return:
        """
        if self.attachment_sketch in integral_flanges or \
                        self.attachment_sketch in optional_flanges:
            return self.F / self.h_o
        elif self.attachment_sketch in loose_flanges:
            return self.F_L / self.h_o
        else:
            raise ValueError("Invalid value %r for property 'Attachment Sketch'" % self.attachment_sketch)

    @property
    def L(self):
        return ((self.t * self.e + 1) / self.T) + ((self.t ** 3) / self.d)

    @property
    def f(self):
        if self.attachment_sketch in loose_flanges or self.hub is None:
            return 1
        else:
            table = Table_2_7_1(self)
            value = table.C_36 / (1 + table.A)
            if value < 1:
                return 1
            else:
                return value

    @property
    def f_actual(self):
        if self.attachment_sketch in loose_flanges or self.hub is None:
            return 1
        else:
            table = Table_2_7_1(self)
            return table.C_36 / (1 + table.A)

    @property
    def S_H_seating(self):
        """
        Longitudinal hub stress
        :return:
        """
        if self.hub is not None:
            return self.f * self.M_o_seating / self.L * (self.hub.g_1 ** 2) * self.B
        else:
            return 0

    @property
    def S_H_operating(self):
        """
        Longitudinal hub stress
        :return:
        """
        if self.hub is not None:
            return self.f * self.M_o_operating / self.L * (self.hub.g_1 ** 2) * self.B
        else:
            return 0

    @property
    def S_R_operating(self):

        if self.hub is not None:
            return ((1.33 * self.t * self.e + 1) * self.M_o_operating) / (self.L * (self.t ** 2) * self.B)
        else:
            return 0

    @property
    def S_R_seating(self):
        if self.hub is not None:
            return ((1.33 * self.t * self.e + 1) * self.M_o_seating) / (self.L * (self.t ** 2) * self.B)
        else:
            return 0

    @property
    def S_T_operating(self):
        if self.hub is not None:
            return ((self.Y * self.M_o_operating) / ((self.t ** 2) * self.B)) - self.Z * self.S_R_operating
        else:
            return (self.Y * self.M_o_operating) / ((self.t ** 2) * self.B)

    @property
    def S_T_seating(self):
        if self.hub is not None:
            return ((self.Y * self.M_o_seating) / ((self.t ** 2) * self.B)) - self.Z * self.S_R_seating
        else:
            return (self.Y * self.M_o_seating) / ((self.t ** 2) * self.B)

    @property
    def K(self):
        return self.A / self.B

    @property
    def T(self):
        num = (self.K ** 2) * (1 + 8.55246 * math.log10(self.K)) - 1
        den = (1.04720 + 1.9448 * (self.K ** 2)) * (self.K - 1)
        return num / den

    @property
    def U(self):
        num = (self.K ** 2) * (1 + 8.55246 * math.log10(self.K)) - 1
        den = 1.36136 * ((self.K ** 2) - 1) * (self.K - 1)
        return num / den

    @property
    def Y(self):
        return (1 / (self.K - 1)) * (0.66845 + 5.71690 * (((self.K ** 2) * math.log10(self.K)) / ((self.K ** 2) - 1)))

    @property
    def Z(self):
        return ((self.K ** 2) + 1) / ((self.K ** 2) - 1)

    @property
    def K_rigidity(self):
        if self.__custom_rigidity_factor is not None and float(self.__custom_rigidity_factor):
            return self.__custom_rigidity_factor
        elif self.attachment_sketch in integral_flanges or \
                        self.attachment_sketch in optional_flanges:
            return 0.3
        elif self.attachment_sketch in loose_flanges:
            return 0.2
        else:
            raise ValueError("Improper value for attribute 'K'")

    @K_rigidity.setter
    def K_rigidity(self, custom_rigidity_factor: float):
        self.__custom_rigidity_factor = float(custom_rigidity_factor)

    @property
    def J_operating(self):
        if self.attachment_sketch in integral_flanges or \
                        self.attachment_sketch in optional_flanges:
            return (52.14 * self.V * self.M_o_operating) / (
                self.L * self.E_f * (self.hub.g_o ** 2) * self.K_rigidity * self.h_o)
        elif self.attachment_sketch in loose_flanges:
            if self.hub is None:
                return (52.14 * self.V_L * self.M_o_operating) / (
                    self.L * self.E_f * (self.hub.g_o ** 2) * self.K_rigidity * self.h_o)
            else:
                return 109.4 * self.M_o_operating / (
                    self.E_f * (self.t ** 3) * self.K_rigidity * math.log(self.K, math.e))
        else:
            raise ValueError("Invalid value %r for attribute 'attachment sketch'" % self.attachment_sketch)

    @property
    def J_seating(self):
        if self.attachment_sketch in integral_flanges or \
                        self.attachment_sketch in optional_flanges:
            return (52.14 * self.V * self.M_o_seating) / (
                self.L * self.E_a * (self.hub.g_o ** 2) * self.K_rigidity * self.h_o)
        elif self.attachment_sketch in loose_flanges:
            if self.hub is None:
                return (52.14 * self.V_L * self.M_o_seating) / (
                    self.L * self.E_a * (self.hub.g_o ** 2) * self.K_rigidity * self.h_o)
            else:
                return 109.4 * self.M_o_seating / (
                    self.E_a * (self.t ** 3) * self.K_rigidity * math.log(self.K, math.e))
        else:
            raise ValueError("Invalid value %r for attribute 'attachment sketch'" % self.attachment_sketch)

    @property
    def PWHT_thickness(self):
        return min(self.t, ((self.A - self.B) * 0.5))

    @property
    def RT_thickness(self):
        return self.PWHT_thickness

    @property
    def PWHT_required(self):
        return self.PWHT_thickness > 3


class Table_2_7_1:
    def __init__(self, parent: Optional[Appendix2FlangeCalcs]=None, g_o=None, g_1=None, h=None, h_o=None):
        self.parent = parent
        if self.parent is None:
            self.g_o = g_o
            self.g_1 = g_1
            self.h = h
            self.h_o = h_o
        else:
            self.g_o = self.parent.hub.g_o
            self.g_1 = self.parent.hub.g_1
            self.h = self.parent.hub.h
            self.h_o = self.parent.hub.h_o

        if not (self.g_o and self.g_1 and self.h and self.h_o):
            raise ValueError("Invalid values passed to table 2-7.1 class.")

    @property
    def A(self):
        return (self.g_1 / self.g_o) - 1

    @property
    def C(self):
        return 43.68 * (self.h / self.h_o) ** 4.0

    @property
    def C_1(self):
        return (1.0 / 3.0) + (self.A / 12.0)

    @property
    def C_2(self):
        return (5.0 / 42.0) + (17.0 * self.A / 336.0)

    @property
    def C_3(self):
        return (1.0 / 210.0) + (self.A / 360.0)

    @property
    def C_4(self):
        return (11.0 / 360.0) + (59.0 * self.A / 5040.0) + ((1.0 + 3.0 * self.A) / self.C)

    @property
    def C_5(self):
        return (1.0 / 90.0) + (5.0 * self.A / 1008.0) - (((1.0 + self.A) ** 3.0) / self.C)

    @property
    def C_6(self):
        return (1 / 120) + (17 * self.A / 5040) + (1 / self.C)

    @property
    def C_7(self):
        return (215 / 2772) + (51 * self.A / 1232) + (60 / 7 + 225 * self.A / 14 + 75 *
                                                      (self.A ** 2) / 7 + 5 * (self.A ** 3) / 2) / self.C

    @property
    def C_8(self):
        return 31 / 6930 + 128 * self.A / 45045 + \
               (6 / 7 + 15 * self.A / 7 + 12 * (self.A ** 2) / 7 + 5 * (self.A ** 3) / 11) / self.C

    @property
    def C_9(self):
        return 533 / 30240 + 653 * self.A / 73920 + \
               (1 / 2 + 33 * self.A / 14 + 39 * (self.A ** 2) / 28 + 25 * (self.A ** 3) / 84) / self.C

    @property
    def C_10(self):
        return 29 / 3780 + 3 * self.A / 704 - \
               (1 / 2 + 33 * self.A / 14 + 81 * (self.A ** 2) / 28 + 13 * (self.A ** 3) / 12) / self.C

    @property
    def C_11(self):
        return 31 / 6048 + 1763 * self.A / 665280 + \
               (1 / 2 + 6 * self.A / 7 + 15 * (self.A ** 2) / 28 + 5 * (self.A ** 3) / 42) / self.C

    @property
    def C_12(self):
        return 1 / 2925 + 71 * self.A / 300300 + \
               (8 / 35 + 18 * self.A / 35 + 156 * (self.A ** 2) / 385 + 6 * (self.A ** 3) / 55) / self.C

    @property
    def C_13(self):
        return 761 / 831600 + 937 * self.A / 1663200 + \
               (1 / 35 + 6 * self.A / 35 + 11 * (self.A ** 2) / 70 + 3 * (self.A ** 3) / 70) / self.C

    @property
    def C_14(self):
        return 197 / 415800 + 103 * self.A / 332640 - \
               (1 / 35 + 6 * self.A / 35 + 17 * (self.A ** 2) / 70 + (self.A ** 3) / 10) / self.C

    @property
    def C_15(self):
        return 233 / 831600 + 97 * self.A / 554400 + \
               (1 / 35 + 3 * self.A / 35 + (self.A ** 2) / 14 + 2 * (self.A ** 3) / 105) / self.C

    @property
    def C_16(self):
        return self.C_1 * self.C_7 * self.C_12 + self.C_2 * self.C_8 * self.C_3 + self.C_3 * self.C_8 * self.C_2 \
               - ((self.C_3 ** 2) * self.C_7 + (self.C_8 ** 2) * self.C_1 + (self.C_2 ** 2) * self.C_12)

    @property
    def C_17(self):
        return (self.C_4 * self.C_7 * self.C_12 + self.C_2 * self.C_8 * self.C_13 + self.C_3 * self.C_8 * self.C_9 -
                (self.C_13 * self.C_7 * self.C_3 +
                 (self.C_8 ** 2) * self.C_4 + self.C_12 * self.C_2 * self.C_9)) / self.C_16

    @property
    def C_18(self):
        return (self.C_5 * self.C_7 * self.C_12 + self.C_2 * self.C_8 * self.C_14 + self.C_3 * self.C_8 * self.C_10 -
                (self.C_14 * self.C_7 * self.C_3 +
                 (self.C_8 ** 2) * self.C_5 + self.C_12 * self.C_2 * self.C_10)) / self.C_16

    @property
    def C_19(self):
        return (self.C_6 * self.C_7 * self.C_12 + self.C_2 * self.C_8 * self.C_15 + self.C_3 * self.C_8 * self.C_11 -
                (self.C_15 * self.C_7 * self.C_3 +
                 (self.C_8 ** 2) * self.C_6 + self.C_12 * self.C_2 * self.C_11)) / self.C_16

    @property
    def C_20(self):
        return (self.C_1 * self.C_9 * self.C_12 + self.C_4 * self.C_8 * self.C_3 + self.C_3 * self.C_13 * self.C_2 -
                (self.C_13 * self.C_8 * self.C_1 +
                 (self.C_3 ** 2) * self.C_9 + self.C_12 * self.C_4 * self.C_2)) / self.C_16

    @property
    def C_21(self):
        return (self.C_1 * self.C_10 * self.C_12 + self.C_5 * self.C_8 * self.C_3 + self.C_3 * self.C_4 * self.C_2 -
                (self.C_14 * self.C_8 * self.C_1 +
                 (self.C_3 ** 2) * self.C_10 + self.C_12 * self.C_5 * self.C_2)) / self.C_16

    @property
    def C_22(self):
        return (self.C_1 * self.C_11 * self.C_12 + self.C_6 * self.C_8 * self.C_3 + self.C_3 * self.C_15 * self.C_2 -
                (self.C_15 * self.C_8 * self.C_1 +
                 (self.C_3 ** 2) * self.C_11 + self.C_12 * self.C_6 * self.C_2)) / self.C_16

    @property
    def C_23(self):
        return (self.C_1 * self.C_7 * self.C_13 + self.C_2 * self.C_9 * self.C_3 + self.C_4 * self.C_8 * self.C_2 -
                (self.C_3 * self.C_7 * self.C_4 +
                 (self.C_2 ** 2) * self.C_13 + self.C_8 * self.C_9 * self.C_1)) / self.C_16

    @property
    def C_24(self):
        return (self.C_1 * self.C_7 * self.C_14 + self.C_2 * self.C_10 * self.C_3 + self.C_5 * self.C_8 * self.C_2 -
                (self.C_3 * self.C_7 * self.C_6 +
                 (self.C_2 ** 2) * self.C_14 + self.C_8 * self.C_10 * self.C_1)) / self.C_16

    @property
    def C_25(self):
        return (self.C_1 * self.C_7 * self.C_15 + self.C_2 * self.C_11 * self.C_3 + self.C_6 * self.C_8 * self.C_2 -
                (self.C_3 * self.C_7 * self.C_6 +
                 (self.C_2 ** 2) * self.C_15 + self.C_8 * self.C_11 * self.C_1)) / self.C_16

    @property
    def C_26(self):
        return (-1) * (self.C / 4) ** 0.25

    @property
    def C_27(self):
        return self.C_20 - self.C_17 - 5 / 12 + self.C_17 * self.C_26

    @property
    def C_28(self):
        return self.C_22 - self.C_19 - (1 / 12) + self.C_19 * self.C_26

    @property
    def C_29(self):
        return -1 * math.sqrt((self.C / 4.0))

    @property
    def C_30(self):
        return -(self.C / 4) ** 0.75

    @property
    def C_31(self):
        return 3 * self.A / 2 - self.C_17 * self.C_30

    @property
    def C_32(self):
        return 0.5 - self.C_19 * self.C_30

    @property
    def C_33(self):
        return 0.5 * self.C_26 * self.C_32 + self.C_28 * self.C_31 * self.C_29 - \
               (0.5 * self.C_30 * self.C_28 + self.C_32 * self.C_27 * self.C_29)

    @property
    def C_34(self):
        return (1 / 12) + self.C_18 - self.C_21 - self.C_18 * self.C_26

    @property
    def C_35(self):
        return - self.C_18 * (self.C / 4) ** 0.75

    @property
    def C_36(self):
        return (self.C_28 * self.C_35 * self.C_29 - self.C_32 * self.C_34 * self.C_29) / self.C_33

    @property
    def C_37(self):
        return (0.5 * self.C_26 * self.C_35 + self.C_34 * self.C_31 * self.C_29 -
                (0.5 * self.C_30 * self.C_34 + self.C_35 * self.C_27 * self.C_29)) / self.C_33

    @property
    def E_1(self):
        return self.C_17 * self.C_36 + self.C_18 + self.C_19 * self.C_37

    @property
    def E_2(self):
        return self.C_20 * self.C_36 + self.C_21 + self.C_22 * self.C_37

    @property
    def E_3(self):
        return self.C_23 * self.C_36 + self.C_24 + self.C_25 * self.C_37

    @property
    def E_4(self):
        return 0.25 + self.C_37 / 12 + self.C_36 / 4 - self.E_3 / 5 - 3 * self.E_2 / 2 - self.E_1

    @property
    def E_5(self):
        return self.E_1 * (0.5 + self.A / 6) + self.E_2 * (0.25 + 11 * self.A / 84) + self.E_3 * \
                                                                                      ((1 / 70) + self.A / 105)

    @property
    def E_6(self):
        return self.E_5 - self.C_36 * ((7 / 120) + self.A / 36 + 3 * self.A / self.C) - \
               (1 / 40) - self.A / 72 - self.C_37 * ((1 / 60) + self.A / 120 + 1 / self.C)

    def input_dict(self):
        return {
            "g_o": self.g_o,
            "g_1": self.g_1,
            "h": self.h,
            "h_o": self.h_o
        }

    def __repr__(self):
        return {
                "A": self.A,
                "C": self.C,
                "C_1": self.C_1,
                "C_2": self.C_2,
                "C_3": self.C_3,
                "C_4": self.C_4,
                "C_5": self.C_5,
                "C_6": self.C_6,
                "C_7": self.C_7,
                "C_8": self.C_8,
                "C_9": self.C_9,
                "E_6": self.E_6,
                "C_10": self.C_10,
                "C_11": self.C_11,
                "C_12": self.C_12,
                "C_13": self.C_13,
                "C_14": self.C_14,
                "C_15": self.C_15,
                "C_16": self.C_16,
                "C_17": self.C_17,
                "C_18": self.C_18,
                "C_19": self.C_19,
                "C_20": self.C_20,
                "C_21": self.C_21,
                "C_22": self.C_22,
                "C_23": self.C_23,
                "C_24": self.C_24,
                "C_25": self.C_25,
                "C_26": self.C_26,
                "C_27": self.C_27,
                "C_28": self.C_28,
                "C_29": self.C_29,
                "C_30": self.C_30,
                "C_31": self.C_31,
                "C_32": self.C_32,
                "C_33": self.C_33,
                "C_34": self.C_34,
                "C_35": self.C_35,
                "C_36": self.C_36,
                "C_37": self.C_37,
                "E_1": self.E_1,
                "E_2": self.E_2,
                "E_3": self.E_3,
                "E_4": self.E_4,
                "E_5": self.E_5
        }


class Units(enum.Enum):
    Imperial = 1
    MKS = 2
    SI = 3


class FacingType(enum.Enum):
    no_facing = 1
    raised_face = 2
    confined_face = 3


class HubGeometry:
    def __init__(self, parent, smallendthickness, largeendthickness, length):
        self.parent = parent
        self.g_o = smallendthickness
        self.g_1 = largeendthickness
        self.h = length


class Table_2_5_2_Sketch(enum.Enum):
    s_1a = 1
    s_1b = 2
    s_1c = 3
    s_1d = 4
    s_2 = 5
    s_3 = 6
    s_4 = 7
    s_5 = 8
    s_6 = 9


class Figure2_4(enum.Enum):
    sketch_1 = 1
    sketch_1a = 2
    sketch_2 = 3
    sketch_2a = 24
    sketch_3 = 4
    sketch_3a = 5
    sketch_4 = 6
    sketch_4a = 7
    sketch_4b = 8
    sketch_4c = 9
    sketch_5 = 10
    sketch_6 = 11
    sketch_6a = 12
    sketch_6b = 13
    sketch_7 = 14
    sketch_8 = 15
    sketch_8a = 16
    sketch_9 = 17
    sketch_9a = 18
    sketch_10 = 19
    sketch_10a = 20
    sketch_11 = 21
    sketch_12 = 22
    sketch_12a = 23


integral_flanges = [
    Figure2_4.sketch_1,
    Figure2_4.sketch_1a,
    Figure2_4.sketch_2,
    Figure2_4.sketch_2a,
    Figure2_4.sketch_3,
    Figure2_4.sketch_3a,
    Figure2_4.sketch_4,
    Figure2_4.sketch_4a,
    Figure2_4.sketch_4b,
    Figure2_4.sketch_4c,
]

loose_flanges = [
    Figure2_4.sketch_5,
    Figure2_4.sketch_6,
    Figure2_4.sketch_6a,
    Figure2_4.sketch_6b,
    Figure2_4.sketch_7
]

optional_flanges = [
    Figure2_4.sketch_8,
    Figure2_4.sketch_8a,
    Figure2_4.sketch_9,
    Figure2_4.sketch_9a,
    Figure2_4.sketch_10,
    Figure2_4.sketch_10a,
    Figure2_4.sketch_11,
    Figure2_4.sketch_12,
    Figure2_4.sketch_12a
]

if __name__ == "__main__":
    go = 5.0
    g1 = 20.0
    length = 32.0
    factor = 63.25  # 800 mm bore
    table = Table_2_7_1(None, go, g1, length, factor)

    print("A is: %s" % table.A)
    print("C is: %s" % table.C)
    print("C_1 is: %s" % table.C_1)
    print("C_2 is: %s" % table.C_2)
    print("C_3 is: %s" % table.C_3)
    print("C_4 is: %s" % table.C_4)
    print("C_5 is: %s" % table.C_5)
    print("C_6 is: %s" % table.C_6)
    print("C_7 is: %s" % table.C_7)
    print("C_8 is: %s" % table.C_8)
    print("C_9 is: %s" % table.C_9)
    print("E_6 is: %s" % table.E_6)
    print("C_10 is: %s" % table.C_10)
    print("C_11 is: %s" % table.C_11)
    print("C_12 is: %s" % table.C_12)
    print("C_13 is: %s" % table.C_13)
    print("C_14 is: %s" % table.C_14)
    print("C_15 is: %s" % table.C_15)
    print("C_16 is: %s" % table.C_16)
    print("C_17 is: %s" % table.C_17)
    print("C_18 is: %s" % table.C_18)
    print("C_19 is: %s" % table.C_19)
    print("C_20 is: %s" % table.C_20)
    print("C_21 is: %s" % table.C_21)
    print("C_22 is: %s" % table.C_22)
    print("C_23 is: %s" % table.C_23)
    print("C_24 is: %s" % table.C_24)
    print("C_25 is: %s" % table.C_25)
    print("C_26 is: %s" % table.C_26)
    print("C_27 is: %s" % table.C_27)
    print("C_28 is: %s" % table.C_28)
    print("C_29 is: %s" % table.C_29)
    print("C_30 is: %s" % table.C_30)
    print("C_31 is: %s" % table.C_31)
    print("C_32 is: %s" % table.C_32)
    print("C_33 is: %s" % table.C_33)
    print("C_34 is: %s" % table.C_34)
    print("C_35 is: %s" % table.C_35)
    print("C_36 is: %s" % table.C_36)
    print("C_37 is: %s" % table.C_37)
    print("E_1 is: %s" % table.E_1)
    print("E_2 is: %s" % table.E_2)
    print("E_3 is: %s" % table.E_3)
    print("E_4 is: %s" % table.E_4)
    print("E_5 is: %s" % table.E_5)
