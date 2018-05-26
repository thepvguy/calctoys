import enum
import math
from Flange.common.Div1Common import Table_2_7_1


class Material:
    def __init__(self, table):
        self.table = table
        self.__stress_data = {}
        self.__get_data()

    def __get_data(self):
        if self.table in [2, 3]:
            self.__stress_data.update({90: 20000})
            self.__stress_data.update({1000: 20000})
        elif self.table in [1]:
            self.__stress_data.update({90: 23000})
            self.__stress_data.update({1000: 23000})

    def allowable(self, temperature):
        if temperature > 90:
            return self.__stress_data[90]
        else:
            return self.__stress_data[1000]


class Gasket:
    def __init__(self, outer_diameter, inner_diameter, gasket_factor, seating_stress):
        self.od = outer_diameter
        self.id = inner_diameter
        self.m = gasket_factor
        self.y = seating_stress


class Bolting:
    def __init__(self, material, diameter, root_area):
        self.material = material
        self.D = diameter
        self.A_r = root_area


class Hub:
    class Sketch(enum.Enum):
        a = enum.auto()
        b = enum.auto()
        c = enum.auto()
        d = enum.auto()

    def __init__(
            self,
            sketch,
            material: Material,
            outside_diameter,
            inner_diameter,
            hub_cross_section_corner_radius,
            small_end_thickness,
            length_of_small_end,
            taper_length,
            neck_length,
            neck_thickness_at_shoulder,
            shoulder_thickness,
            shoulder_height,
            transition_angle,
            shoulder_angle,
            friction_angle
    ):
        self.sketch = sketch
        self.material = material
        self.A = outside_diameter
        self.B = inner_diameter
        self.r = hub_cross_section_corner_radius
        self.g_0 = small_end_thickness
        self.h_small = length_of_small_end
        self.h = taper_length
        self.h_n = neck_length
        self.g_1 = neck_thickness_at_shoulder
        self.T = shoulder_thickness
        self.g_2 = shoulder_height
        self.alpha = transition_angle
        self.phi = shoulder_angle
        self.mu = friction_angle


    @property
    def f(self):
        table = Table_2_7_1(g_o=self.g_0, g_1=self.g_1, h=self.h, h_o=(self.B * self.g_0) ** 0.5)
        return max(table.C_36 / (1 + table.A), 1.0)

    @property
    def h_2(self):
        return self.T - 0.5 * (self.g_2 * math.tan(math.radians(self.phi)))

    @property
    def g_bar(self):
        num = self.T * (self.g_1 ** 2) + self.h_2 * self.g_2 * (2 * self.g_1 + self.g_2)
        den = (2 * (self.T * self.g_1 + self.h_2 * self.g_2))
        return num / den

    @property
    def h_bar(self):
        num = (self.T ** 2) * self.g_1 + self.g_2 * self.h_2 ** 2
        den = 2 * (self.T * self.g_1 + self.h_2 * self.g_2)
        return num / den

    @property
    def I_h(self):
        term1 = (1 / 3) * self.g_1 * self.T ** 3
        term2 = (1 / 3) * self.g_2 * self.h_2 ** 3
        term3 = (self.g_2 * self.h_2 + self.g_1 * self.T) * self.h_bar ** 2
        return term1 + term2 - term3

    @property
    def N(self):
        return self.B + 2 * self.g_1

    @property
    def M_o_factor(self):
        term1 = 1.818 / ((self.B * self.g_1) ** 0.5)
        term2 = (3.305 * self.I_h) / ((0.5 * self.B + self.g_bar) * self.g_1 ** 2)
        return 1 + term1 * (self.T - self.h_bar + term2)

    def Z(self, is_operating):
        if is_operating:
            return self.phi - self.mu
        else:
            return self.phi + self.mu

    def get_result(self):
        result = {
            "f": self.f,
            "h_2": self.h_2,
            "g_bar": self.g_bar,
            "h_bar": self.h_bar,
            "I_h": self.I_h,
            "N": self.N
        }
        return result


class Clamp:
    class Sketch(enum.Enum):
        a = enum.auto()
        b = enum.auto()
        c = enum.auto()

    def __init__(
            self,
            sketch,
            material: Material,
            bolt_circle_radius,
            clamp_inside_diameter,
            clamp_width,
            effective_clamp_thickness,
            effective_clamp_gap,
            corner_radius,
            distance_from_bolt_circle_to_clamp_od,
            effective_lip_length,
            lug_height,
            lug_width
    ):
        self.sketch = sketch
        self.material = material
        self.B_c = bolt_circle_radius
        self.C_i = clamp_inside_diameter
        self.C_w = clamp_width
        self.C_t = effective_clamp_thickness
        self.C_g = effective_clamp_gap
        self.r = corner_radius
        self.L_a = distance_from_bolt_circle_to_clamp_od
        self.l_c = effective_lip_length
        self.L_h = lug_height
        self.L_w = lug_width

    @property
    def A_1(self):
        return (self.C_w - 2 * self.C_t) * self.C_t

    @property
    def A_2(self):
        return 1.571 * self.C_t ** 2

    @property
    def A_3(self):
        return (self.C_w - self.C_g) * self.l_c

    @property
    def A_c(self):
        return self.A_1 + self.A_2 + self.A_3

    @property
    def X(self):
        return (((0.5 * self.C_w - (1 / 3) * self.C_t) * (self.C_t ** 2)) - 0.5 * (self.C_w - self.C_g) * self.l_c ** 2) / self.A_c

    @property
    def e_b(self):
        return self.B_c - 0.5 * self.C_i - self.l_c - self.X

    @property
    def I_c(self):
        return ((1 / 3) * self.A_1 + 0.25 * self.A_2)*(self.C_t ** 2) + (1 / 3)*(self.A_3 * self.l_c ** 2) - self.A_c * self.X ** 2

    def get_result(self):
        result = {
            "A_1": self.A_1,
            "A_2": self.A_2,
            "A_3": self.A_3,
            "A_c": self.A_c,
            "X": self.X,
            "e_b": self.e_b,
            "I_c": self.I_c
        }
        return result


class DesignCondition:
    def __init__(
            self,
            is_operating,
            temperature,
            ambient_temperature,
            pressure,
            bolting_is_controlled,
            has_retainer
    ):
        self.is_operating = is_operating
        self.temperature = temperature
        self.ambient_temperature = ambient_temperature
        self.pressure = pressure
        self.bolting_controlled = bolting_is_controlled
        self.has_retainer = has_retainer


class HubAndClamp:
    def __init__(self, design_data: DesignCondition, bolting: Bolting, gasket: Gasket, hub: Hub, clamp: Clamp):
        self.design_data = design_data
        self.bolts = bolting
        self.gasket = gasket
        self.hub = hub
        self.clamp = clamp

        self.do_sanity_checks()

    def do_sanity_checks(self):
        if not (self.clamp.sketch in [Clamp.Sketch.a, Clamp.Sketch.b, Clamp.Sketch.c]):
            ValueError("Invalid value for clamp sketch")

    # -- BEGIN BOLT LOADS --

    @property
    def N_bolts(self):
        if self.clamp.sketch == Clamp.Sketch.b:
            return 2
        else:
            return 4

    @property
    def G(self):
        # TODO -- make this more detailed. Implement table 2-5.2
        return 0.5 * (self.gasket.od - self.gasket.id) + self.gasket.id

    @property
    def b(self):
        # TODO -- make this more accurate. Implement table 2-5.2
        return 0.5 * (self.gasket.od - self.gasket.id)

    @property
    def H(self):
        return 0.785 * self.design_data.pressure * self.G ** 2

    @property
    def H_p(self):
        return 2 * self.b * math.pi * self.G * self.gasket.m * self.design_data.pressure

    @property
    def W_m_1(self):
        return 0.637 * (self.H + self.H_p) * math.tan(math.radians(self.hub.phi - self.hub.mu))

    @property
    def H_m(self):
        #  TODO -- implement custom input for total axial seating load
        return math.pi * self.b * self.G * self.gasket.y

    @property
    def W_m_2(self):
        return 0.637 * self.H_m * math.tan(math.radians(self.hub.phi + self.hub.mu))

    @property
    def W_m_3(self):
        return 0.637 * (self.H + self.H_p) * math.tan(math.radians(self.hub.phi + self.hub.mu))

    @property
    def S_a(self):
        return self.bolts.material.allowable(self.design_data.ambient_temperature)

    @property
    def S_b(self):
        return self.bolts.material.allowable(self.design_data.temperature)

    @property
    def A_m_1(self):
        return self.W_m_1 / (2 * self.S_b)

    @property
    def A_m_2(self):
        return self.W_m_2 / (2 * self.S_a)

    @property
    def A_m_3(self):
        return self.W_m_3 / (2 * self.S_a)

    @property
    def A_m_L(self):
        return max(self.A_m_1, self.A_m_2, self.A_m_3)

    @property
    def A_b_L(self):
        return self.bolts.A_r * self.N_bolts / 2

    @property
    def W(self):
        if self.design_data.is_operating:
            return self.W_m_1
        else:
            if self.design_data.bolting_controlled:
                return 2 * self.A_m_L * self.bolts.material.allowable(self.design_data.ambient_temperature)
            else:
                return (self.A_m_L + self.A_b_L) * self.bolts.material.allowable(self.design_data.ambient_temperature)

    # -- END BOLT LOADS --
    # -- START HUB MOMENTS --

    @property
    def C(self):
        return (self.hub.A + self.clamp.C_i) * 0.5

    @property
    def H_D(self):
        return 0.785 * self.design_data.pressure * self.hub.B ** 2

    @property
    def h_D(self):
        return 0.5 * (self.C - (self.hub.B + self.hub.g_1))

    @property
    def M_D(self):
        return self.H_D * self.h_D

    @property
    def M_F(self):
        return self.H_D * 0.5 * (self.hub.g_1 - self.hub.g_0)

    @property
    def H_G(self):
        return ((1.571 * self.W) / (math.tan(math.radians(self.hub.phi + self.hub.mu)))) - (self.H + self.H_p)

    @property
    def h_G(self):
        # TODO: dig into this and figure out what else it should be
        return 0

    @property
    def H_T(self):
        return self.H - self.H_D

    @property
    def h_T(self):
        return 0.5 * (self.C - 0.5 * (self.G + self.hub.B))

    @property
    def M_G(self):
        return self.H_G * self.h_G

    @property
    def M_P(self):
        return math.pi * self.design_data.pressure * self.hub.B * self.hub.T * (0.5 * self.hub.T - self.hub.h_bar)

    @property
    def M_R(self):
        return 1.571 * self.W * (self.hub.h_bar - self.hub.T + 0.5 * (self.C - self.hub.N) * math.tan(math.radians(self.hub.phi)))

    @property
    def M_T(self):
        return self.H_T * self.h_T

    @property
    def M_o(self):
        if self.design_data.is_operating:
            return self.M_D + self.M_G + self.M_T + self.M_F + self.M_P + self.M_R
        else:
            return (0.785 * self.W * (self.C - self.G)) / (math.tan(math.radians(self.hub.phi + self.hub.mu)))

    @property
    def M_H(self):
        return self.M_o / self.hub.M_o_factor

    @property
    def Q(self):
        return 1.818 * self.M_H / ((self.hub.B * self.hub.g_1) ** 0.5)

    # -- END HUB MOMENTS --
    # -- START HUB STRESSES --

    def __Z(self):
        return self.hub.Z(self.design_data.is_operating)

    @property
    def S_1(self):
        term1 = ((self.design_data.pressure * self.hub.B ** 2)/(4 * self.hub.g_1 * (self.hub.B + self.hub.g_1)))
        term2 = (1.91 * self.M_H) / ((self.hub.g_1 ** 2) * (self.hub.B + self.hub.g_1))
        return self.hub.f * (term1 + term2)

    @property
    def S_2(self):
        return self.design_data.pressure * (((self.hub.N ** 2) + (self.hub.B ** 2)) / ((self.hub.N ** 2) - (self.hub.B ** 2)))

    @property
    def S_3(self):
        return (0.75 * self.W) / (self.hub.T * (self.hub.B + self.hub.g_1 * 2) * math.tan(math.radians(self.__Z())))

    @property
    def S_4(self):
        return (0.477 * self.Q) / (self.hub.g_1 * (self.hub.B + self.hub.g_1))

    # -- END HUB STRESSES --
    # -- START CLAMP STRESSES --

    @property
    def l_m(self):
        return self.clamp.l_c - (self.C - self.clamp.C_i) * 0.5

    @property
    def S_5(self):
        term1 = self.W / (2 * self.C * math.tan(math.radians(self.__Z())))
        term2 = 1 / self.clamp.C_t
        term3 = (3 * (self.clamp.C_t + 2 * self.l_m)) / self.clamp.C_t ** 2
        return term1 * (term2 + term3)

    @property
    def S_6(self):
        term1 = self.W / 2
        term2 = (1 / self.clamp.A_c)
        term3 = (abs(self.clamp.e_b) * (self.clamp.C_t - self.clamp.X)) / self.clamp.I_c
        return term1 * (term2 + term3)

    @property
    def S_7(self):
        return (1.5 * self.W) / ((self.clamp.C_w - self.clamp.C_g) * self.C * math.tan(math.radians(self.__Z())))

    @property
    def S_8(self):
        return (3 * self.W * self.clamp.L_a) / (self.clamp.L_w * self.clamp.L_h ** 2)

    @property
    def S_9(self):
        return self.W / ((self.hub.A - self.clamp.C_i) * self.C * math.tan(math.radians(self.__Z())))

    # -- END CLAMP STRESSES --
    # -- START ACCEPTIBILITY CRITERIA --

    @property
    def S_H(self):
        if self.design_data.is_operating:
            return self.hub.material.allowable(self.design_data.temperature)
        else:
            return self.hub.material.allowable(self.design_data.ambient_temperature)

    @property
    def S_C(self):
        if self.design_data.is_operating:
            return self.clamp.material.allowable(self.design_data.temperature)
        else:
            return self.clamp.material.allowable(self.design_data.ambient_temperature)

    @property
    def max_S1(self):
        return 1.5 * self.S_H
        
    @property
    def max_S2(self):
        return self.S_H
    
    @property
    def max_S3(self):
        return 0.8 * self.S_H
        
    @property
    def max_S4(self):
        return 0.8 * self.S_H
        
    @property
    def max_S5(self):
        return 1.5 * self.S_C
        
    @property
    def max_S6(self):
        return 1.5 * self.S_C
        
    @property
    def max_S7(self):
        return self.S_C
        
    @property
    def max_S8(self):
        return self.S_C
        
    @property
    def max_S9(self):
        return 1.6 * min(self.S_H, self.S_C)
    
    @property
    def S1_is_acceptable(self):
        return self.S_1 <= self.max_S1
        
    @property
    def S2_is_acceptable(self):
        return self.S_2 <= self.max_S2
        
    @property
    def S3_is_acceptable(self):
        return self.S_3 <= self.max_S3
        
    @property
    def S4_is_acceptable(self):
        return self.S_4 <= self.max_S4
        
    @property
    def S5_is_acceptable(self):
        return self.S_5 <= self.max_S5
    
    @property
    def S6_is_acceptable(self):
        return self.S_6 <= self.max_S6
        
    @property
    def S7_is_acceptable(self):
        return self.S_7 <= self.max_S7
        
    @property
    def S8_is_acceptable(self):
        return self.S_8 <= self.max_S8
        
    @property
    def S9_is_acceptable(self):
        return self.S_9 <= self.max_S9

    def get_result(self):
        result = {
            "N_bolts": self.N_bolts,
            "G": self.G,
            "b": self.b,
            "H": self.H,
            "H_p": self.H_p,
            "W_m_1": self.W_m_1,
            "H_m": self.H_m,
            "W_m_2": self.W_m_2,
            "W_m_3": self.W_m_3,
            "S_a": self.S_a,
            "S_b": self.S_b,
            "A_m_1": self.A_m_1,
            "A_m_2": self.A_m_2,
            "A_m_3": self.A_m_3,
            "A_m_L": self.A_m_L,
            "A_b_L": self.A_b_L,
            "W": self.W,
            "C": self.C,
            "H_D": self.H_D,
            "h_D": self.h_D,
            "M_D": self.M_D,
            "M_F": self.M_F,
            "H_G": self.H_G,
            "h_G": self.h_G,
            "H_T": self.H_T,
            "h_T": self.h_T,
            "M_G": self.M_G,
            "M_P": self.M_P,
            "M_R": self.M_R,
            "M_T": self.M_T,
            "M_o": self.M_o,
            "M_H": self.M_H,
            "Q": self.Q,
            "Z": self.__Z(),
            "hub": self.hub.get_result(),
            "clamp": self.clamp.get_result(),
            "S_1": self.S_1,
            "S_2": self.S_2,
            "S_3": self.S_3,
            "S_4": self.S_4,
            "l_m": self.l_m,
            "S_5": self.S_5,
            "S_6": self.S_6,
            "S_7": self.S_7,
            "S_8": self.S_8,
            "S_9": self.S_9,
            "S_H": self.S_H,
            "S_C": self.S_C,
            "max_S1": self.max_S1,
            "max_S2": self.max_S2,
            "max_S3": self.max_S3,
            "max_S4": self.max_S4,
            "max_S5": self.max_S5,
            "max_S6": self.max_S6,
            "max_S7": self.max_S7,
            "max_S8": self.max_S8,
            "max_S9": self.max_S9,
            "S1_is_acceptable": self.S1_is_acceptable,
            "S2_is_acceptable": self.S2_is_acceptable,
            "S3_is_acceptable": self.S3_is_acceptable,
            "S4_is_acceptable": self.S4_is_acceptable,
            "S5_is_acceptable": self.S5_is_acceptable,
            "S6_is_acceptable": self.S6_is_acceptable,
            "S7_is_acceptable": self.S7_is_acceptable,
            "S8_is_acceptable": self.S8_is_acceptable,
            "S9_is_acceptable": self.S9_is_acceptable
        }
        return result


class Assembly:
    def __init__(
            self,
            design_data,
            clamp,
            left_hub,
            right_hub,
            bolting
    ):
        # TODO : Evaluate assembly and operating conditions
        # TODO : Evaluate different right and left clamps, maybe
        self.design_data = design_data
        self.clamp = clamp
        self.left_hub = left_hub
        self.right_hub = right_hub
        self.bolting = bolting
