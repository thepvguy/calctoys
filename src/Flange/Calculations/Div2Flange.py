import enum
import math


class Units(enum.Enum):
    customary = enum.auto()
    metric = enum.auto()

class GasketFacingSketch(enum.Enum):
    sketch_1a = enum.auto()
    sketch_1b = enum.auto()
    sketch_1c = enum.auto()
    sketch_1d = enum.auto()
    sketch_2 = enum.auto()
    sketch_3 = enum.auto()
    sketch_4 = enum.auto()
    sketch_5 = enum.auto()
    sketch_6 = enum.auto()


class GasketFacingSketchColumn(enum.Enum):
    column_1 = enum.auto()
    column_2 = enum.auto()


class FlangeType(enum.Enum):
    integral = enum.auto()
    loose = enum.auto()
    reverse_integral = enum.auto()
    reverse_loose = enum.auto()


class IntegralType(enum.Enum):
    # Figure 4.16.1
    no_hub = enum.auto()
    straight_hub = enum.auto()
    # Figure 4.16.2
    hub_type_1 = enum.auto()
    hub_type_2 = enum.auto()
    hub_type_3 = enum.auto()
    # Figure 4.16.3
    nut_stop_detail_A = enum.auto()
    nut_stop_detail_B = enum.auto()
    nut_stop_detail_C = enum.auto()
    nut_stop_detail_D = enum.auto()
    # Figure 4.16.4
    nut_stop_large = enum.auto()


class LooseType(enum.Enum):
    # Figure 4.16.5
    with_hub = enum.auto()
    without_hub = enum.auto()
    # Figure 4.16.6
    lap_with_hub = enum.auto()
    lap_no_hub = enum.auto()


#  GASKETS
class Div2FlangeGasketParameters:
    def __init__(
            self,
            units,
            facing_sketch,
            sketch_column,
            inner_diameter,
            outer_diameter,
            factor_m,
            factor_y,
            nubbin_width=-1,
            thickness=-1

    ):
        self.units = units
        self.facing_sketch = facing_sketch
        self.sketch_column = sketch_column
        self.inner_diameter = inner_diameter
        self.outer_diameter = outer_diameter
        self.factor_m = factor_m
        self.factor_y = factor_y
        self.nubbin_width = nubbin_width
        self.thickness = thickness


class Div2FlangeGasket:
    def __init__(self, params: Div2FlangeGasketParameters):
        self.p = params

    @property
    def outer_diameter(self):
        return self.p.outer_diameter

    @property
    def m(self):
        return self.p.factor_m

    @property
    def y(self):
        return self.p.factor_y

    def N(self):
        return (self.p.outer_diameter - self.p.inner_diameter) / 2

    def C_ul(self):
        if self.p.units == Units.customary:
            return 1
        elif self.p.units == Units.metric:
            return 25.4
        else:
            raise ValueError("Invalid value passed to units type!")

    def __b_0_lower_bound(self):
        if self.p.units == Units.customary:
            return 0.25
        elif self.p.units == Units.metric:
            return 6
        else:
            raise ValueError("Invalid value passed to units type!")

    def b(self):
        if self.b_0() <= self.__b_0_lower_bound():
            return self.b_0()
        else:
            return 0.5 * self.C_ul() * (self.b_0() / self.C_ul()) ** 0.5

    def b_0(self):
        if self.p.facing_sketch == GasketFacingSketch.sketch_1a or self.p.facing_sketch == GasketFacingSketch.sketch_1b:
            return self.N() / 2
        elif self.p.facing_sketch == GasketFacingSketch.sketch_1c or self.p.facing_sketch == GasketFacingSketch.sketch_1d:
            return min((self.p.nubbin_width + self.p.thickness) / 2, (self.p.nubbin_width + self.N()) / 4)
        elif self.p.facing_sketch == GasketFacingSketch.sketch_2:
            if self.p.facing_sketch == GasketFacingSketchColumn.column_1:
                return (self.p.nubbin_width + self.N()) / 4
            else:
                return (self.p.nubbin_width + 3 * self.N()) / 8
        elif self.p.facing_sketch == GasketFacingSketch.sketch_3 or self.p.facing_sketch == GasketFacingSketch.sketch_5:
            if self.p.facing_sketch == GasketFacingSketchColumn.column_1:
                return self.N() / 4
            else:
                return 3 * self.N() / 8
        elif self.p.facing_sketch == GasketFacingSketch.sketch_4:
            if self.p.facing_sketch == GasketFacingSketchColumn.column_1:
                return 3 * self.N() / 8
            else:
                return 7 * self.N() / 16
        elif self.p.facing_sketch == GasketFacingSketch.sketch_6:
            if self.p.facing_sketch == GasketFacingSketchColumn.column_1:
                return self.p.nubbin_width / 8
        else:
            raise ValueError("No valid value of b_0 for sketch %s and column %s")


#  BOLTS
class Div2Bolt:
    def __init__(self, nominal_diameter, root_area, allowable_stress_operating, allowable_stress_ambient):
        self.nominal_diameter = nominal_diameter
        self.root_area = root_area
        self.allowable_stress_operating = allowable_stress_operating
        self.allowable_stress_ambient = allowable_stress_ambient



#  FLANGE
class Div2FlangeDesignConditions:
    def __init__(
            self,
            internal_pressure,
            allowable_stress,
            axial_load,
            bending_moment,
            units,
            flange_type,
            integral_type,
            loose_type
    ):

        self.internal_pressure = internal_pressure
        self.allowable_stress = allowable_stress
        self.axial_load = axial_load
        self.bending_moment = bending_moment
        self.units = units
        self.flange_type = flange_type
        self.integral_type = integral_type
        self.loose_type = loose_type


class Div2FlangeHubGeometry:
    def __init__(
            self,
            length,
            small_end_thickness,
            large_end_thickness
    ):
        self.h = length
        self.g_0 = small_end_thickness
        self.g1 = large_end_thickness


class Div2FlangeGeometry:
    def __init__(
            self,
            outer_diameter,
            bore_diameter,
            bolt_circle_diameter,
            thickness,
            number_of_bolts,
            facing_inner_diameter,
            facing_outer_diameter,
            facing_is_confined,
            number_of_splits = 0,
            hub_geometry=None
    ):
        self.A = outer_diameter
        self.B = bore_diameter
        self.C = bolt_circle_diameter
        self.t = thickness
        self.num_bolts = number_of_bolts
        self.facing_inner_diameter = facing_inner_diameter
        self.facing_outer_diameter = facing_outer_diameter
        self.confined_face = facing_is_confined
        self.fs_num = number_of_splits
        self.hub_geometry = hub_geometry


class Div2Flange:
    def __init__(
            self,
            conditions: Div2FlangeDesignConditions,
            geometry: Div2FlangeGeometry,
            bolts: Div2Bolt,
            gasket: Div2FlangeGasket):
        self.conditions = conditions
        self.geometry = geometry
        self.bolts = bolts
        self.gasket = gasket

    def G_c(self):
        return min(self.gasket.outer_diameter, self.geometry.facing_outer_diameter)

    def G(self):
        return self.G_c() - self.gasket.b()

    def W_o(self):
        load = 0.785 * (self.G() ** 2) * self.conditions.internal_pressure
        if self.gasket.m != 0 and self.gasket.y != 0:  # Non self energizing
            load += self.gasket.b() * math.pi * self.G() * self.gasket.m * self.conditions.internal_pressure

        return load

    def W_g(self):
        return 0.5 * (self.A_m() + self.A_b()) * self.bolts.allowable_stress_ambient

    def F_A(self):
        if self.conditions.axial_load < 0:
            return abs(self.conditions.internal_pressure)
        else:
            return self.conditions.axial_load

    def M_E(self):
        return abs(self.conditions.bending_moment)

    def C_us(self):
        if self.conditions.units == Units.customary:
            return 1
        elif self.conditions.units == Units.metric:
            return 6.894757 * 10 ** -3

    def W_gs(self):
        if self.gasket.m == 0 and self.gasket.y == 0:
            return 0
        else:
            return math.pi * self.gasket.b() * self.G() * self.C_us() * self.gasket.y

    def __A_m_ambient(self):
        return self.W_gs() / self.bolts.allowable_stress_ambient

    def __A_m_operating(self):
        return (self.W_o() + self.F_A() + (4 * self.M_E() / self.G())) / self.bolts.allowable_stress_operating

    def ambient_bolt_load_governs(self):
        return self.__A_m_ambient() > self.__A_m_operating()

    def A_m(self):
        if self.ambient_bolt_load_governs():
            return self.__A_m_ambient()
        else:
            return self.__A_m_operating()

    def A_b(self):
        return self.bolts.root_area * self.geometry.num_bolts

    def A_b_ok(self):
        return self.A_b() >= self.A_m()

    def H_D(self):
        return 0.785 * (self.geometry.B ** 2) * self.conditions.internal_pressure

    def h_D(self):
        if self.conditions.flange_type == FlangeType.integral:
            return (self.geometry.C - self.geometry.B - self.geometry.hub_geometry.g_1) / 2
        elif self.is_loose_flange():
            return 0.5 * (self.geometry.C - self.geometry.B)
        elif self.conditions.flange_type == FlangeType.reverse_integral:
            return 0.5 * (self.geometry.C + self.geometry.hub_geometry.g_1 - 2 * self.geometry.hub_geometry.g_o - self.geometry.B)
        else:
            raise ValueError("Invalid configuration given for h_d")

    def H(self):
        return 0.785 * (self.G() ** 2) * self.conditions.internal_pressure

    def H_T(self):
        return self.H() - self.H_D()

    def h_T(self):
        t = self.conditions.flange_type

        if t == FlangeType.integral:
            return 0.5 * (0.5 * (self.geometry.C - self.geometry.B) + self.h_G())
        elif t == FlangeType.loose:
            if self.conditions.loose_type in [LooseType.lap_no_hub, LooseType.lap_with_hub]:
                return (self.geometry.C - self.G()) * 0.5
            elif self.conditions.loose_type in [LooseType.with_hub, LooseType.without_hub]:
                return (self.h_D() + self.h_G()) * 0.5
            else:
                raise ValueError("Invalid value for loose type configuration given for h_T")
        elif self.is_reverse_flange():
            return 0.5 * (self.geometry.C - 0.5 * (self.geometry.B + self.G()))

        else:
            raise ValueError("Invalid flange type given for value h_T")

    def H_G(self):
        return self.W_o() - self.H()

    def h_G(self):
        return 0.5 * (self.geometry.C - self.G())

    def B_star(self):
        RuntimeWarning("B_star has not been looked at, please implement this; results using this value should be carefully examined")
        return self.geometry.B

    def is_loose_flange(self):
        if self.conditions.flange_type in [FlangeType.loose, FlangeType.reverse_loose]:
            return True
        elif self.conditions.flange_type in [FlangeType.integral, FlangeType.reverse_integral]:
            return False
        else:
            raise ValueError("Invalid flange type")

    def is_reverse_flange(self):
        if self.conditions.flange_type in [FlangeType.reverse_integral, FlangeType.reverse_loose]:
            return True
        elif self.conditions.flange_type in [FlangeType.integral, FlangeType.loose]:
            return False
        else:
            raise ValueError("Invalid flange type")

    def has_hub(self):
        if self.geometry.hub_geometry is not None:
            return False
        else:
            return True

    def __hub_factors_applicable(self, is_reverse):
        if self.conditions.flange_type == FlangeType.integral or (self.is_loose_flange() and self.has_hub()):
            return not is_reverse
        elif self.conditions.flange_type == FlangeType.reverse_integral or self.conditions.loose_type in [LooseType.with_hub, LooseType.lap_with_hub]:
            return is_reverse
        else:
            raise ValueError("Invalid configuration used when querying hub factor applicability")

    def K(self):
        if self.__hub_factors_applicable(False):
            return self.geometry.A / self.geometry.B
        elif self.__hub_factors_applicable(True):
            return self.geometry.A / self.B_star()
        else:
            raise ValueError("Invalid configuration")

    def Y(self):
        return (1 / (self.K() - 1)) * (0.66845 + 5.71690 * ((self.K()**2) * math.log10(self.K())) / ((self.K() ** 2) - 1))

    def T(self):
        num = (self.K() ** 2) * (1 + 8.55246 * math.log10(self.K())) - 1
        den = (1.04720 + 1.9448 * self.K() ** 2) * (self.K() - 1)
        return num / den

    def U(self):
        num = (self.K() ** 2) * (1 + 8.55246 * math.log10(self.K())) - 1
        den = 1.36136 * ((self.K() ** 2) - 1) * (self.K() - 1)
        return num / den

    def Z(self):
        return ((self.K() ** 2) + 1) / ((self.K() ** 2) - 1)

    def h_o(self):
        return math.sqrt(self.geometry.B * self.geometry.hub_geometry.g_o)

    def e(self):
        if not self.is_loose_flange():
            return self.F() / self.h_o()
        else:
            if self.has_hub():
                return self.F_L() / self.h_o()
            else:
                raise ValueError("Invalid flange configuraation for e")

    def d(self):
        num = self.U() * (self.geometry.hub_geometry.g_o ** 2) * self.h_o()

        if not self.is_loose_flange():
            den = self.V()
        else:
            if self.has_hub():
                den = self.V_L()
            else:
                raise ValueError("Invalid flange configuraation for d")

        return num / den

    def L(self):
        return ((self.geometry.t * self.e() + 1) / self.T()) + ((self.geometry.t ** 3) / (self.d()))

    def X_g(self):
        return self.geometry.hub_geometry.g_1 / self.geometry.hub_geometry.g_0

    def X_h(self):
        if self.is_reverse_flange():
            return self.geometry.hub_geometry.h / self.h_or()
        else:
            return self.geometry.hub_geometry.h / self.h_o()

    def alpha_r(self):
        return (1 / self.K() ** 2) * (1 + ((0.668 * (self.K() + 1)) / (self.Y())))

    def Y_r(self):
        return self.Y() * self.alpha_r()

    def T_r(self):
        return ((self.Z() + 0.3) / (self.Z() - 0.3)) * self.alpha_r() * self.T()

    def U_r(self):
        return self.alpha_r() * self.U()

    def h_or(self):
        return math.sqrt(self.geometry.A * self.geometry.hub_geometry.g_o)

    def e_r(self):
        if not self.is_loose_flange():
            return self.F() / self.h_or()
        else:
            if self.has_hub():
                return self.F_L() / self.h_or()
            else:
                raise ValueError("Invalid flange configuraation for e")

    def d_r(self):
        num = self.U_r() * (self.geometry.hub_geometry.g_o ** 2) * self.h_or()

        if not self.is_loose_flange():
            den = self.V()
        else:
            if self.has_hub():
                den = self.V_L()
            else:
                raise ValueError("Invalid flange configuraation for d")

        return num / den

    def L_r(self):
        return ((self.geometry.t * self.e_r() + 1) / (self.T_r())) + ((self.geometry.t ** 3) / self.d_r())

    def F(self):
        ln_xh = math.log(self.X_h(), math.e)
        ln_xg = math.log(self.X_g(), math.e)
        ret = 0.897697
        ret += 9.5257 * (10 ** -3) * ln_xh
        ret += 0.123586 * ln_xg ** 2
        ret += 0.0358580 * ln_xh ** 2
        ret -= 0.194422 * ln_xg * ln_xh
        ret -= 0.0181259 * ln_xg ** 3
        ret += 0.0129360 * ln_xh ** 3
        ret -= 0.0377693 * ln_xg * ln_xh ** 2
        ret += 0.0273791 * ln_xh * ln_xg ** 2

        return ret

    def F_L(self):
        lnxg = math.log(self.X_g(), math.e)
        lnxh = math.log(self.X_h(), math.e)

        num = 0.941074
        num += 0.176139 * lnxg
        num -= 0.188556 * lnxh
        num += 0.0689847 * lnxg ** 2
        num += 0.523798 * lnxg ** 2
        num -= 0.513894 * lnxg * lnxh

        den = 1
        den += 0.379392 * lnxg
        den += 0.184520 * lnxh
        den -= 0.00605208 * lnxg ** 2
        den -= 0.00358934 * lnxh ** 2
        den += 0.110179 * lnxg * lnxh

        return num / den

    def V(self):
        if self.X_h() < 0.1 or self.X_h() > 2:
            raise ValueError("X_h is out of range for V")
        elif 0.1 <= self.X_h() <= 0.5:
            ret = 0.500244
            ret += (0.227914 / self.X_g())
            ret -= 1.87071 * self.X_h()
            ret -= (0.344410 / (self.X_g() ** 2))
            ret += 2.49189 * self.X_h() ** 2
            ret += 0.873446 * (self.X_h() / self.X_g())
            ret += 0.189953 / (self.X_g() ** 3)
            ret -= 1.06082 * self.X_h() ** 3
            ret -= 1.49970 * ((self.X_h() ** 2) / self.X_g())
            ret += 0.719413 * (self.X_h() / (self.X_g() ** 2))

        elif 0.5 < self.X_h() <= 2:
            xg_inv = 1 / self.X_g()
            xh_inv = 1/ self.X_h()
            ret = 0.0144868
            ret -= 0.135977 * xg_inv
            ret -= 0.0461919 * xh_inv
            ret += 0.560718 * xg_inv ** 2
            ret += 0.0529829 * xh_inv ** 2
            ret += 0.244313 * xg_inv * xh_inv
            ret += 0.113929 * xg_inv ** 3
            ret -= 0.00928265 * xh_inv ** 3
            ret -= 0.0266293 * xg_inv * xh_inv ** 2
            ret -= 0.217008 * xh_inv * xg_inv ** 2
        else:
            raise ValueError("Invalid value for X_h")

        return ret

    def V_L(self):
        xh = self.X_h()

        if xh < 0.1 or xh > 2.0:
            raise ValueError("X_h is out of range")

        xg = self.X_g()
        lnxh = math.log(xh, math.e)
        lnxg = math.log(xg, math.e)
        xh_inv = 1 / xh
        xg_inv = 1 / xg

        if 0.1 <= xh <= 0.25:
            ret = 6.57683
            ret -= 0.115516 * xg
            ret += 1.39499 * (xg ** 0.5) * lnxg
            ret += 0.307340 * lnxg ** 2
            ret -= 8.30849 * xg ** 0.5
            ret += 2.62307 * lnxg
            ret += 0.239498 * xh * lnxh
            ret -= 2.96125 * lnxh
            ret += (7.035052 * 10 ** -4) * xh_inv

            ret = math.e ** ret

        elif 0.25 < xh <= 0.5:
            ret = 1.56323
            ret -= 1.80696 * lnxg
            ret -= 1.33458 * xh_inv
            ret += 0.276415 * lnxg ** 2
            ret += 0.417135 * xh_inv ** 2
            ret += 1.39511 * lnxg * xh_inv
            ret += 0.0137129 * lnxg ** 3
            ret += 0.0943597 * xh_inv ** 3
            ret -= 0.402096 * lnxg * xh_inv ** 2
            ret -= 0.101619 * xh_inv * lnxg ** 2

        elif 0.5 < xh <= 1.0:

            ret = -0.0213643
            ret -= 0.0763597 * xg_inv
            ret += 0.1029900 * xh_inv
            ret += 0.725776 * xg_inv ** 2
            ret -= 0.160603 * xh_inv ** 2
            ret -= 0.0918061 * xg_inv * xh_inv
            ret += 0.472277 * xg_inv ** 3
            ret += 0.0873530 * xh_inv ** 3
            ret += 0.527487 * xg_inv * xh_inv ** 2
            ret -= 0.980209 * xh_inv * xg_inv ** 2

        elif 1.0 < xh <= 2.0:
            ret = 7.96687 * (10 ** -3)
            ret -= 0.220518 * xg_inv
            ret += 0.0602652 * xh_inv
            ret += 0.619818 * xg_inv ** 2
            ret -= 0.223212 * xh_inv ** 2
            ret += 0.421920 * xg_inv * xh_inv
            ret += 0.0950195 * xg_inv ** 3
            ret += 0.209813 * xh_inv ** 3
            ret -= 0.158821 * xg_inv * xh_inv ** 2
            ret -= 0.242056 * xh_inv * xg_inv ** 2
        else:
            raise ValueError("Invalid value for X_h")

        return ret

    def f(self):
        if self.is_loose_flange():
            return 1

        num = 0.0927779
        num -= 0.0336633 * self.X_g()
        num += 0.964176 * self.X_g() ** 2
        num += 0.0566286 * self.X_h()
        num += 0.347074 * self.X_h() ** 2
        num -= 4.18699 * self.X_h() ** 3

        den = 1
        den -= 5.96093 * (10 ** -3) * self.X_g()
        den += 1.62904 * self.X_h()
        den += 3.49329 * self.X_h() ** 2
        den += 1.39052 * self.X_h() ** 3

        return max(1.0, num / den)

    def I(self):
        if self.conditions.flange_type == FlangeType.integral:
            if self.conditions.integral_type == IntegralType.no_hub:
                raise ValueError("Table 4.16.7 does not provide calculations for integral flanges with no hub")
            elif self.conditions.integral_type in []:
                return (0.0874 * self.L() * (self.geometry.hub_geometry.g_o ** 2) * self.h_o() * self.geometry.B) / self.V()
            else:
                raise ValueError("Invalid value set for flange integral type in I")
        elif self.conditions.flange_type == FlangeType.loose:
            if self.conditions.loose_type in [LooseType.with_hub, LooseType.lap_with_hub]:
                return (0.0874 * self.L() * (self.geometry.hub_geometry.g_o ** 2) * self.h_o() * self.geometry.B) / self.V_L()
            elif self.conditions.loose_type in [LooseType.without_hub, LooseType.lap_no_hub]:
                return (1/6) * (self.geometry.B * (self.geometry.t ** 3) * math.log(self.K(), math.e))
            else:
                raise ValueError("Invalid value for loose type in value I()")
        else:
            raise ValueError("Invalid value for flange type in I()")

    def A_R(self):
        return 0.5 * (self.geometry.A - self.geometry.B)

    def G_avg(self):
        return 0.5 * (self.geometry.hub_geometry.g_0 + self.geometry.hub_geometry.g_1)

    def A_A(self):
        if self.geometry.t >= self.G_avg():
            return self.A_R()
        else:
            return self.geometry.hub_geometry.h + self.geometry.t

    def B_B(self):
        if self.geometry.t >= self.G_avg():
            return self.geometry.t
        else:
            return self.G_avg()

    def C_C(self):
        if self.geometry.t >= self.G_avg():
            return self.geometry.hub_geometry.h
        else:
            return self.A_R() - self.G_avg()

    def D_DG(self):
        if self.geometry.t >= self.G_avg():
            return self.G_avg()
        else:
            return self.geometry.t

    def K_AB(self):
        return (self.A_A() * self.B_B() ** 3) * ((1/3) - 0.21 * (self.B_B() / self.A_A()) * (1 - (1/12) * (self.B_B() / self.A_A()) ** 4))

    def K_CD(self):
        return (self.C_C() * self.D_DG() ** 3) * ((1/3) - 0.105 * (self.D_DG()/self.C_C()) * (1 - (1 / 192) * (self.D_DG() / self.C_C()) ** 4))

    def I_p(self):
        if self.has_hub():
            return self.K_AB() + self.K_CD()
        else:
            return self.A_R() * self.geometry.t ** 3 * ((1/3) - 0.21 * (self.geometry.t / self.A_R()) * (1 - (1/12) * (self.geometry.t / self.A_R()) ** 4))

    def M_oe(self):
        t1 = 4 * self.M_E()
        t2 = self.I() / (0.3846 * self.I_p() + self.I())
        t3 = self.h_D() / (self.geometry.C - 2 * self.h_D())

        return t1 * t2 * t3 + self.F_A() * self.h_D()

    def F_S(self):
        if self.geometry.fs_num == 0:
            return 1.0
        elif self.geometry.fs_num == 1:
            return 2.0
        elif self.geometry.fs_num == 2:
            return 0.75
        else:
            raise ValueError("invalid number of splits specified for flange")

    def B_s(self):
        return (0.25 * math.pi * self.geometry.B ** 2) / self.geometry.num_bolts

    def B_sc(self):
        return max(1.0, math.sqrt(self.B_s() / (2 * self.bolts.nominal_diameter + self.geometry.t)))

    def B_s_max(self):
        return 2 * self.bolts.nominal_diameter + ((6 * self.geometry.t) / (self.gasket.m + 0.5))

    def __M_o_int(self):
        return abs(((self.H_D() * self.h_D() + self.H_T() * self.h_T() + self.H_G() * self.h_G()) * self.B_sc() + self.M_oe()) * self.F_S())

    def __M_o_ext(self):
        return abs((self.H_D() * (self.h_D() - self.h_G()) + self.H_T() * (self.h_T() - self.h_G()) + self.M_oe()) * self.F_S())

    def M_o(self):
        if self.conditions.internal_pressure >=0:
            return self.__M_o_int()
        else:
            return self.__M_o_ext()

    def __M_g_ext(self):
        return self.W_g() * self.h_G() * self.F_S()

    def __M_g_int(self):
        return 0.5 * (self.W_g() * (self.geometry.C - self.G()) * self.B_sc() * self.F_S())

    def M_g(self):
        if self.conditions.internal_pressure >=0:
            return self.__M_g_int()
        else:
            return self.__M_g_ext()

    def S_H_operating(self):
        if self.is_loose_flange() and not self.has_hub():
            return None

        if self.is_reverse_flange():
            return self.f() * self.M_o() / (self.L_r() * self.B_star() * self.geometry.hub_geometry.g_1 ** 2)
        else:
            return self.f() * self.M_o() / (self.L() * self.B() * self.geometry.hub_geometry.g_1 ** 2)

    def S_H_seating(self):
        if self.is_loose_flange() and not self.has_hub():
            return None

        if self.is_reverse_flange():
            return self.f() * self.M_g() / (self.L_r() * self.B_star() * self.geometry.hub_geometry.g_1 ** 2)
        else:
            return self.f() * self.M_g() / (self.L() * self.geometry.B * self.geometry.hub_geometry.g_1 ** 2)

    def S_R_operating(self):
        if self.is_loose_flange() and not self.has_hub():
            return None

        if self.is_reverse_flange():
            return (1.33 * self.geometry.t * self.e_r() + 1) * self.M_o() / (self.L_r() * self.B_star() * self.geometry.t ** 2)
        else:
            return (1.33 * self.geometry.t * self.e() + 1) * self.M_o() / (self.L() * self.geometry.B * self.geometry.t ** 2)

    def S_R_seating(self):
        if self.is_loose_flange() and not self.has_hub():
            return None

        if self.is_reverse_flange():
            return (1.33 * self.geometry.t * self.e_r() + 1) * self.M_g() / (self.L_r() * self.B_star() * self.geometry.t ** 2)
        else:
            return (1.33 * self.geometry.t * self.e() + 1) * self.M_g() / (self.L() * self.geometry.B * self.geometry.t ** 2)

    def S_T_operating(self):

        if self.is_reverse_flange():
            if self.is_loose_flange() and not self.has_hub():
                return self.Y() * self.M_o() / (self.B_star() * self.geometry.t ** 2)
        else:
            if (self.is_loose_flange() and self.has_hub()) or not self.is_loose_flange():
                return (self.Y() * self.M_o() / (self.geometry.B * self.geometry.t ** 2)) - self.Z() * self.S_R_operating()
            elif self.is_loose_flange() and not self.has_hub():
                return self.Y() * self.M_o() / (self.geometry.B * self.geometry.t ** 2)
            else:
                raise ValueError("What configuration do you have?")

    def S_T_seating(self):

        if self.is_reverse_flange():
            if self.is_loose_flange() and not self.has_hub():
                return self.Y() * self.M_g() / (self.B_star() * self.geometry.t ** 2)
        else:
            if (self.is_loose_flange() and self.has_hub()) or not self.is_loose_flange():
                return (self.Y() * self.M_g() / (self.geometry.B * self.geometry.t ** 2)) - self.Z() * self.S_R_seating()
            elif self.is_loose_flange() and not self.has_hub():
                return self.Y() * self.M_g() / (self.geometry.B * self.geometry.t ** 2)
            else:
                raise ValueError("What configuration do you have?")

    def S_T_1_operating(self):
        if self.is_reverse_flange() and self.has_hub():
            t1 = (self.Y_r() * self.M_o()) / (self.B_star() * self.geometry.t ** 2)
            t2 = (self.Z() * self.S_R_operating() * (0.67 * self.geometry.t * self.e_r() + 1)) / (1.33 * self.geometry.t * self.e_r() + 1)
            return t1 - t2

    def S_T_1_seating(self):
        if self.is_reverse_flange() and self.has_hub():
            t1 = (self.Y_r() * self.M_g()) / (self.B_star() * self.geometry.t ** 2)
            t2 = (self.Z() * self.S_R_seating() * (0.67 * self.geometry.t * self.e_r() + 1)) / (1.33 * self.geometry.t * self.e_r() + 1)
            return t1 - t2

    def S_T_2_operating(self):
        if self.is_reverse_flange() and self.has_hub():
            num = ((0.67 * self.geometry.t * self.e_r() + 1) * 2 * self.K() ** 2)
            den = (((self.K() ** 2) - 1) * self.L_r())
            term1 = (self.Y() - (num / den))
            term2 = (self.M_o() / self.B_star() * self.geometry.t ** 2)
            return term1 * term2

    def S_T_2_seating(self):
        if self.is_reverse_flange() and self.has_hub():
            num = ((0.67 * self.geometry.t * self.e_r() + 1) * 2 * self.K() ** 2)
            den = (((self.K() ** 2) - 1) * self.L_r())
            term1 = (self.Y() - (num / den))
            term2 = (self.M_g() / self.B_star() * self.geometry.t ** 2)
            return term1 * term2

    def check_S_H_operating(self):
        pass



    def getResult(self):
        return {
            "bolt_loads": {
                "b": self.gasket.b(),
                "G": self.G(),
                "W_0": self.W_o(),
                "W_gs": self.W_gs(),
                "ambient_bolt_load_governs": self.ambient_bolt_load_governs(),
                "A_m": self.A_m(),
                "A_b": self.A_b(),
                "A_b_ok": self.A_b_ok()
            },

            "factors": {

            },

            "forces": {
                "H_D": self.H_D(),
                "H": self.H(),
                "H_T": self.H_T(),
                "H_G": self.H_G()
            },

            "moments": {

            },

            "stresses": {

            },

            "acceptance": {

            },

            "rigidity": {

            }

        }



