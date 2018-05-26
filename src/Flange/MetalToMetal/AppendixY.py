import math
import enum
import abc
import Flange.common.Div1Common as common

class Material:
    def __init__(self, table):
        self.table=table


class Bolt:
    def __init__(
            self,
            material,
            nominal_diameter,
            root_area,
            is_stud
    ):
        self.material = material
        self.nominal_diameter = nominal_diameter
        self.root_area = root_area
        self.is_stud = is_stud

    def E(self, temperature):
        return 23000


class Gasket:
    def __init__(
            self,
            inner_diameter,
            outer_diameter,
            thickness,
            gasket_factor,
            seating_stress
    ):
        self.id = inner_diameter
        self.od = outer_diameter
        self.t = thickness
        self.m = gasket_factor
        self.y = seating_stress

    @property
    def G(self):
        return 0.5 * (self.od - self.id) + self.id


class Hub:
    def __init__(
            self,
            small_end_thickness,
            large_end_thickness,
            overall_length
    ):
        self.g_o = small_end_thickness
        self.g_1 = large_end_thickness
        self.h = overall_length

class Spacer:
    def __init__(
            self,
            outer_diameter,
            inner_diameter,
            thickness,
            material
    ):
        self.od = outer_diameter
        self.id = inner_diameter
        self.t = thickness
        self.material = material


class AbstractFlange(metaclass=abc.ABCMeta):
    def __init__(
            self,
            design_data,
            material,
            attachment_sketch,
            bolting,
            outer_diameter,
            inner_diameter,
            bolt_circle,
            bolt_hole_diameter,
            number_of_bolts,
            thickness,
            rim_contact_diameter,
            gasket,
            hub=None
    ):
        super().__init__()
        self.design_data = design_data
        self.material = material
        self.bolting = bolting
        self.A = outer_diameter
        self.B = inner_diameter
        self.C = bolt_circle
        self.D = bolt_hole_diameter
        self.n = number_of_bolts
        self.t = thickness
        self.h_C = rim_contact_diameter
        self.gasket = gasket
        self.attachment_sketch = attachment_sketch
        self.hub = hub

    @property
    def category(self):
        if self.attachment_sketch in common.optional_flanges or self.attachment_sketch in common.integral_flanges:
            return FlangeCategory.cat1
        elif self.attachment_sketch in common.loose_flanges and self.hub is not None:
            return FlangeCategory.cat2
        else:
            return FlangeCategory.cat3

    @property
    def AR_bar(self):
        return (self.n * self.D) / (math.pi * self.C)

    @property
    def h_o(self):
        return (self.B * self.hub.g_o) ** 0.5

    @property
    def f(self):
        table = common.Table_2_7_1(
            g_o=self.hub.g_o,
            g_1=self.hub.g_1,
            h=self.hub.h,
            h_o=self.h_o
        )
        return max(1, table.C_36 / (1 + self.A))

    @property
    def B_1(self):
        if self.hub is None:
            return self.B
        elif self.f == 1:
            return self.B + self.hub.g_o
        else:
            return self.B + self.hub.g_1

    @property
    def a(self):
        return (self.A + self.C) / (2 * self.B_1)

    @property
    def r_B(self):
        term1 = 4 / (1 - self.AR_bar ** 2) ** 0.5
        term2 = ((1 + self.AR_bar) / (1 - self.AR_bar)) ** 0.5
        return (1 / self.n) * (term1 * math.atan(term2) - math.pi - 2 * self.AR_bar)

    @property
    def G(self):
        return self.gasket.G

    @property
    def R(self):
        # TODO: take attachment sketch into account
        return 0.5 * (self.C - self.B) - self.hub.g_1

    @property
    def h_D(self):
        table = common.Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_D

    @property
    def h_T(self):
        table = common.Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_T

    @property
    def h_G(self):
        table = common.Table_2_6(
            attachment_sketch=self.attachment_sketch,
            inner_diameter=self.B,
            bolt_circle=self.C,
            large_end_hub_thickness=self.hub.g_1,
            gasket_reaction_diameter=self.G,
            distance_from_bolt_circle_to_hub=self.R
        )
        return table.h_G

    @property
    def H_D(self):
        return 0.785 * self.design_data.P * self.B ** 2

    @property
    def H(self):
        return 0.785 * (self.G ** 2) * self.design_data.P

    @property
    def H_T(self):
        return self.H - self.H_D

    @property
    def beta(self):
        return (self.C + self.B_1) / (2 * self.B_1)

    @property
    def J_S(self):
        return (1 / self.B_1) * (((2 * self.h_D) / self.beta) + (self.h_C / self.a)) + math.pi * self.r_B

    @property
    def J_P(self):
        return (1 / self.B_1) * ((self.h_D / self.beta) + (self.h_C / self.a)) + math.pi * self.r_B

    @property
    def C_1(self):
        num = - (0.784 - 1.567 * self.J_S * math.log10(self.A / self.B_1))
        den = 1 + 1.3 * self.J_S
        return num / den

    @property
    def K(self):
        return self.A / self.B

    @property
    def T(self):
        return common.Figure_2_7_1(self.K).T

    @property
    def U(self):
        return common.Figure_2_7_1(self.K).U

    @property
    def Y(self):
        return common.Figure_2_7_1(self.K).Y

    @property
    def Z(self):
        return common.Figure_2_7_1(self.K).Z

    @property
    def V(self):
        return common.Table_2_7_1(g_o=self.hub.g_o, g_1=self.hub.g_1, h=self.hub.h, h_o=self.h_o).V

    @property
    def V_L(self):
        return common.Table_2_7_1(g_o=self.hub.g_o, g_1=self.hub.g_1, h=self.hub.h, h_o=self.h_o).V_L

    @property
    def F(self):
        return common.Table_2_7_1(g_o=self.hub.g_o, g_1=self.hub.g_1, h=self.hub.h, h_o=self.h_o).F

    @property
    def F_L(self):
        return common.Table_2_7_1(g_o=self.hub.g_o, g_1=self.hub.g_1, h=self.hub.h, h_o=self.h_o).F_L

    @property
    def H_G(self):
        NotImplementedError("")  # TODO
        return 0

    @property
    def E(self):
        return self.material.E(self.design_data.temperature)

    @property
    def M_P(self):
        return self.H_D * self.h_D + self.H_T * self.h_T * self.H_G * self.h_G

    @property
    def r_E(self):
        return self.E / self.bolting.E(self.design_data.temperature)


    # -- BEGIN CLASS SPECIFIC EQUATIONS; POSSIBLE REFACTOR NEEDED --

    @property
    @abc.abstractmethod
    def F_prime(self):
        NotImplementedError("F_prime")
        return 0

    @property
    @abc.abstractmethod
    def l(self):
        NotImplementedError("l")
        return 0

    @property
    @abc.abstractmethod
    def M_S(self):
        NotImplementedError("M_S")
        return 0

    @property
    @abc.abstractmethod
    def E_theta_B(self):
        NotImplementedError("theta_b")
        return 0

    @property
    @abc.abstractmethod
    def H_C(self):
        NotImplementedError("H_C")
        return 0

    @property
    @abc.abstractmethod
    def W_m_1(self):
        NotImplementedError("W_m_1")
        return 0

    @property
    @abc.abstractmethod
    def sigma_b(self):
        NotImplementedError("sigma_b")
        return 0

    @property
    @abc.abstractmethod
    def S_i(self):
        NotImplementedError("S_i")
        return 0

    @property
    @abc.abstractmethod
    def S_R_BC(self):
        NotImplementedError("S_R_BC")
        return 0

    @property
    @abc.abstractmethod
    def S_R(self):
        NotImplementedError("S_R")
        return 0

    @property
    @abc.abstractmethod
    def S_T(self):
        NotImplementedError("S_T")
        return 0

    @property
    @abc.abstractmethod
    def S_H(self):
        NotImplementedError("S_H")
        return 0


class FlangeClass1(AbstractFlange):
    def __init__(
            self,
            design_data,
            material,
            attachment_sketch,
            bolting,
            outer_diameter,
            inner_diameter,
            bolt_circle,
            bolt_hole_diameter,
            number_of_bolts,
            thickness,
            rim_contact_diameter,
            gasket,
            spacer=None,
            hub=None):

        super().__init__(
            design_data,
            material,
            attachment_sketch,
            bolting,
            outer_diameter,
            inner_diameter,
            bolt_circle,
            bolt_hole_diameter,
            number_of_bolts,
            thickness,
            rim_contact_diameter,
            gasket,
            hub
        )

        self.spacer = spacer

    @property
    def t_s(self):
        if self.spacer is not None:
            return self.spacer.thickness
        else:
            return 0

    @property
    def l(self):
        factor = 1 if self.bolting.is_stud else 0.5
        return 2 * self.t + self.t_s + factor * self.bolting.nominal_diameter

    @property
    def M_S(self):
        # Equation 7
        num = self.J_P * self.F_prime * self.M_P
        den = (self.t ** 3) + self.J_S * self.F_prime
        return num / den

    @property
    def E_theta_B(self):
        # Equation 8
        return (5.46 / (math.pi * self.t ** 3)) * (self.J_S * self.M_S + self.J_S * self.M_P)

    @property
    def H_C(self):
        # Equation 9
        return (self.M_P + self.M_S) / self.h_C

    @property
    def W_m_1(self):
        # Equation 10
        return self.H + self.H_G + self.H_C

    @property
    def sigma_b(self):
        # Equation 11
        return self.W_m_1 / (self.bolting.root_area * self.n)

    @property
    def S_i(self):
        num = 1.159 * (self.M_P + self.M_S) * self.h_C ** 2
        den = (self.a * self.l * self.r_E * self.B_1)
        return self.sigma_b - (num / den)

    @property
    def S_R_BC(self):
        # Equation 13
        num = 6 * (self.M_P + self.M_S)
        den = (self.t ** 2) * (math.pi * self.C - self.n * self.D)
        return num / den

    @property
    def S_R(self):
        if self.category == FlangeCategory.cat3:
            return 0

        term2 = self.M_S / (math.pi * self.B_1 * self.t ** 2)

        if self.category == FlangeCategory.cat1:
            term1 = (2 * self.F * self.t) / (self.h_o + self.F * self.t)
        else:  # self.category == FlangeCategory.cat2:
            term1 = (2 * self.F_L * self.t) / (self.h_o + self.F_L * self.t)

        return -(term1 + 6) * term2

    @property
    def S_T(self):
        term1 = self.t * self.E_theta_B / self.B_1

        if self.category == FlangeCategory.cat3:
            return term1

        if self.category == FlangeCategory.cat1:
            term2 = (2 * self.F * self.t * self.Z) / (self.h_o + self.F * self.t)
        else:  # self.category == FlangeCategory.cat2:
            term2 = (2 * self.F_L * self.t * self.Z) / (self.h_o + self.F_L * self.t)

        term3 = self.M_S / (math.pi * self.B_1 * self.t ** 2)

        return term1 + (term2 - 1.8) * term3

    @property
    def S_H(self):
        if self.category == FlangeCategory.cat3:
            return 0

        if self.category == FlangeCategory.cat1:
            v = self.V
        else:  # self.category == FlangeCategory.cat2:
            v = self.V_L

        return (self.h_o * self.E_theta_B * self.f) / (0.91 * ((self.hub.g_1 / self.hub.g_o) ** 2) * self.B_1 * v)


class FlangeCategory(enum.Enum):
    cat1 = enum.auto()
    cat2 = enum.auto()
    cat3 = enum.auto()


class FlangeClass(enum.Enum):
    cls1 = enum.auto()
    cls2 = enum.auto()
    cls3 = enum.auto()