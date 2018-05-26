import math

from Flange.common.Div1Common import Table_2_7_1, Units, FacingType, Table_2_5_2_Sketch, integral_flanges, \
    loose_flanges, optional_flanges, Table_2_6, Figure_2_7_1


class Appendix2Params:
    def __init__(self, other):
        # TODO
        self.internal_pressure = None
        self.bolt_ambient_allowable_stress = None
        self.bolt_operating_allowable_stress = None
        self.flange_operating_allowable_stress = None
        self.nozzle_neck_operating_allowable = None
        self.flange_ambient_modulus_of_elasticity = None
        self.flange_operating_modulus_of_elasticity = None
        self.custom_Wm1_active = None
        self.custom_Wm1 = None
        self.custom_Wm2_active = None
        self.custom_Wm2 = None
        self.gasket_OD = None
        self.gasket_ID = None
        self.gasket_m = None
        self.gasket_y = None
        self.gasket_thickness = None
        self.facing_sketch = None
        self.facing_column = None
        self.w = None
        self.nubbin_height = None
        self.raised_face_type = None
        self.rf_dia = None
        self.flange_OD = None
        self.flange_ID = None
        self.flange_thickness = None
        self.num_bolts = None
        self.bolt_circle_diameter = None
        self.hub = None
        self.bolt_dismeter = None
        self.bolt_root_area = None
        self.attachment_sketch = None
        self.custom_rigidity = None
        self.units = None


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
        self.hub = params.hub

        self.a = params.bolt_dismeter
        self.boltRootArea = params.bolt_root_area

        self.attachment_sketch = params.attachment_sketch
        self.__custom_rigidity_factor = params.custom_rigidity
        self.units = params.units

    @classmethod
    def from_components(cls, design_condition, geometry, materials):
        pass

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

        table = Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_D

    @property
    def h_G(self):
        """
        Radial distance between gasket load reaction and bolt circle
        :return:
        """
        table = Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_G

    @property
    def h_T(self):
        """
        radial distance between bot circle and H_t
        :return:
        """
        table = Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_T

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
        return Table_2_7_1(self).V

    @property
    def V_L(self):
        assert (self.hub is not None)
        return Table_2_7_1(self).V_L

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
        return Table_2_7_1(self).F

    @property
    def F_L(self):
        assert (self.hub is not None)
        return Table_2_7_1(self).F_L

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
            return Table_2_7_1(self).f(True)
        else:
            return Table_2_7_1(self).f(False)

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
        return Figure_2_7_1(self.K).T

    @property
    def U(self):
        return Figure_2_7_1(self.K).U

    @property
    def Y(self):
        return Figure_2_7_1(self.K).Y

    @property
    def Z(self):
        return Figure_2_7_1(self.K).Z

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


class HubGeometry:
    def __init__(self, parent, smallendthickness, largeendthickness, length):
        self.parent = parent
        self.g_o = smallendthickness
        self.g_1 = largeendthickness
        self.h = length


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
