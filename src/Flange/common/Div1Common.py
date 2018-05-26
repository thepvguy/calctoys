import enum
import math

from typing import Optional

from Flange.Traditional.Appendix_2 import Appendix2FlangeCalcs


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
    def F(self):
        if self.g_1 == self.g_o:
            return 0.908920
        else:
            return self.E_6 / (((self.C / 2.73) ** 0.25) * (((1 + self.A) ** 3) / self.C))

    @property
    def F_L(self):
        num1 = self.C_18 * (0.5 + (self.A / 6))
        num2 = self.C_21 * (0.25 + (11 * self.A / 84))
        num3 = self.C_24 * ((1 / 70) + (self.A / 105))
        num4 = (1 / 40) + (self.A / 72)
        den = ((self.C / 2.73) ** 0.25) * (((1 + self.A) ** 3) / self.C)
        return (num1 + num2 + num3 - num4) / den

    @property
    def V(self):
        if self.g_1 == self.g_o:
            return 0.550103
        else:
            return self.E_4 / (((2.73 / self.C) ** 0.25) * ((1 + self.A) ** 3))

    @property
    def V_L(self):
        return (0.25 - 0.2 * self.C_24 - 1.5 * self.C_21 - self.C_18) / (
            ((2.73 / self.C) ** 0.25) * (1 + self.A) ** 3)

    def f(self, is_loose):
        if is_loose:
            return 1
        else:
            return max(1, self.C_36 / (1 + self.A))

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


class Table_2_6:
    def __init__(
            self,
            attachment_sketch,
            inner_diameter,
            bolt_circle,
            large_end_hub_thickness,
            gasket_reaction_diameter,
            distance_from_bolt_circle_to_hub
    ):
        self.attachment_sketch = attachment_sketch
        self.B = inner_diameter
        self.C = bolt_circle
        self.g_1 = large_end_hub_thickness
        self.G = gasket_reaction_diameter
        self.R = distance_from_bolt_circle_to_hub

    @property
    def h_D(self):
        """
        Radial distance between bolt circle and H_d radius
        :return:
        """
        if self.attachment_sketch in integral_flanges:
            return self.R + 0.5 * self.g_1
        elif self.attachment_sketch in loose_flanges or self.attachment_sketch in optional_flanges:
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
            return 0.5 * (self.R + self.g_1 + self.h_G)
        elif self.attachment_sketch in loose_flanges or self.attachment_sketch in optional_flanges:

            if self.attachment_sketch in [Table_2_5_2_Sketch.sketch_1, Table_2_5_2_Sketch.sketch_1a]:
                return 0.5 * (self.C - self.G)
            else:
                return 0.5 * (self.h_D + self.h_G)
        else:
            raise ValueError("Invalid value '%r' for attribute 'attachment_sketch'" % self.attachment_sketch)


class Figure_2_7_1:
    def __init__(self, K):
        self.K = K

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
