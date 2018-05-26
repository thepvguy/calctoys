import math

from _UHX_common import PitchType


class UHX11Params:
    def __init__(
            self,
            radius_to_outermost_tube_hole_center,
            nominal_tube_OD,
            nominal_tube_wall_thickness,
            modulus_of_elasticity_of_tubes_at_tubesheet_design_temperature,
            modulus_of_elasticity_for_tubesheet_material_at_tubesheet_design_temperature,
            allowable_stress_of_tubes_at_tubesheet_design_temperature,
            allowable_stress_for_tubesheet_material_at_tubesheet_design_temperature,
            tube_pitch,
            area_of_untubed_lanes,
            tube_side_pass_partition_groove_depth,
            tubesheet_corrosion_tube_side,
            tubesheet_thickness,
            expanded_depth_of_tube_in_tubesheet,
            pitch_type
    ):
        self.r_o = radius_to_outermost_tube_hole_center
        self.d_t = nominal_tube_OD
        self.t_t = nominal_tube_wall_thickness
        self.Et_T = modulus_of_elasticity_of_tubes_at_tubesheet_design_temperature
        self.E = modulus_of_elasticity_for_tubesheet_material_at_tubesheet_design_temperature
        self.St_T = allowable_stress_of_tubes_at_tubesheet_design_temperature
        self.S = allowable_stress_for_tubesheet_material_at_tubesheet_design_temperature
        self.p = tube_pitch
        self.A_L = area_of_untubed_lanes
        self.h_g = tube_side_pass_partition_groove_depth
        self.c_t = tubesheet_corrosion_tube_side
        self.h = tubesheet_thickness
        self.l_tx = expanded_depth_of_tube_in_tubesheet
        self.pitchType = pitch_type


class UHX11:
    def __init__(self, params):

        self.r_o = params.r_o
        self.d_t = params.d_t
        self.t_t = params.t_t
        self.Et_T = params.Et_T
        self.E = params.E
        self.St_T = params.St_T
        self.S = params.S
        self.p = params.p
        self.A_L = params.A_L
        self.h_g = params.h_g
        self.c_t = params.c_t
        self.h = params.h
        self.l_tx = params.l_tx
        self.pitchType = params.pitchType

    def rho(self):
        """

        :return: Tube expansion depth ratio
        """
        return max(0, min(1, self.l_tx / self.h))


    def D_o(self):
        """

        :return: Equivalent outer tube limit circle
        """
        return 2.0 * self.r_o + self.d_t

    def mu(self):
        """

        :return: basic ligament efficiency for shear
        """
        return (self.p - self.d_t) / self.p

    def d_star(self):
        """

        :return: Effective tube hole diameter
        """
        return max(
            (self.d_t - 2.0 * self.t_t * (self.Et_T / self.E) * (self.St_T / self.S) * self.rho()),
            (self.d_t - 2.0 * self.t_t)
        )

    def p_star(self):
        """

        :return: Effective tube pitch
        """
        return self.p / math.sqrt(1 - ((4.0 / (math.pi * self.D_o() ** 2)) * min(self.A_L, 4.0 * self.D_o() * self.p)))

    def mu_star(self):
        """

        :return: Effective ligament efficiency for bending
        """
        return (self.p_star() - self.d_star()) / self.p_star()

    def hg_prime(self):
        """

        :return: Effective tube side pass partition groove depth
        """
        return max(self.h_g - self.c_t, 0)

    @staticmethod
    def __interpolate(x1, y1, x2, y2, x):
        """

        :param x1: Starting X
        :param y1: Starting Y
        :param x2: Ending X
        :param y2: Ending Y
        :param x: Evaluate at this X
        :return: The interpolated Y value for the given two input points
        """
        return y1 + (x - x1) * ((y2 - y1) / (x2 - x1))

    def __hpRatio(self):
        """

        :return: The ratio of the thickness vs the pitch
        """
        return self.h / self.p

    def __alpha(self, number):
        # TODO: abstract __alpha and __beta further to avoid code duplication
        """

        :param number: The subscript of the corresponding alpha coefficient (e.g. 0 for alpha_0)
        :return: Interpolated alpha of the given subscript
        """
        # alpha is zero indexed, let's keep it that way for external users.
        number += 1

        # This table may change every few years; not often enough to be put in a database,
        # and should definitely be under revision control
        # [h / p, alpha_0, alpha_1, alpha_2, alpha_3, alpha_4]
        alphaTableTriangular = [
            [0.10, 0.0353, 1.2502, -0.0491, 0.3604, -0.6100],
            [0.25, 0.0135, 0.9910, 1.0080, -1.0498, 0.0184],
            [0.50, 0.0054, 0.5279, 3.0461, -4.3657, 1.9435],
            [2.00, -0.0029, 0.2126, 3.9906, -6.1730, 3.4307]
        ]

        alphaTableSquare = [
            [0.10, 0.0676, 1.5756, -1.2119, 1.7715, -1.2628],
            [0.25, 0.0250, 1.9251, -3.5230, 6.9830, -5.0017],
            [0.50, 0.0394, 1.3024, -1.1041, 2.8714, -2.3994],
            [2.00, 0.0372, 1.0314, -0.6402, 2.6201, -2.1929]
        ]
        if self.pitchType in [PitchType.TRIANGLE, PitchType.ROTATED_TRIANGLE]:
            alphaTable = alphaTableTriangular
        elif self.pitchType in [PitchType.SQUARE, PitchType.ROTATED_SQUARE]:
            alphaTable = alphaTableSquare
        else:
            raise ValueError("Invalid pitch type <<%r>> given, which is not in pitch enum PitchType" % self.pitchType)


        alphaOut = None
        hpRatio = self.__hpRatio()

        for i in range(len(alphaTable)):
            alphaRow = alphaTable[i]
            if i == 0:
                if hpRatio <= alphaRow[0]:
                    alphaOut = alphaRow[number]

            elif i == len(alphaTable) and not hpRatio < alphaRow[0]:
                betaOut = betaRow

            else:
                if hpRatio < alphaRow[0]:
                    alphaOut = self.__interpolate(alphaTable[i - 1][0], alphaTable[i - 1][number], alphaTable[i][0],
                                                 alphaTable[i][number], hpRatio)

        return alphaOut

    def __beta(self, number):
        """

        :param number: The subscript of the corresponding beta coefficient (e.g. 0 for beta_0)
        :return: Interpolated beta of the given subscript
        """
        # alpha is zero indexed, let's keep it that way for external users.
        number += 1

        # This table may change every few years; not often enough to be put in a database,
        # and should definitely be under revision control
        # [h / p, alpha_0, alpha_1, alpha_2, alpha_3, alpha_4]
        betaTableTriangular = [
            [0.10, -0.0958, 0.6209, -0.8683, 2.1099, -1.6831],
            [0.15, 0.8897, -9.0855, 36.1435, -59.5425, 35.8223],
            [0.25, 0.7439, -4.4989, 12.5779, -14.2092, 5.7822],
            [0.50, 0.9100, -4.8901, 12.4325, -12.7039, 4.4298],
            [1.00, 0.9923, -4.8759, 12.3572, -13.7214, 5.7629],
            [2.00, 0.9966, -4.1978, 9.0478, -7.9955, 2.2398]
        ]

        betaTableSquare = [
            [0.10, -0.0791, 0.6008, -0.3468, 0.4858, -0.3606],
            [0.15, 0.3345, -2.8420, 10.9709, -15.8994, 8.3516],
            [0.25, 0.4296, -2.6350, 8.6864, -11.5227, 5.8544],
            [0.50, 0.3636, -0.8057, 2.0463, -2.2902, 1.1862],
            [1.00, 0.3527, -0.2842, 0.4354, -0.901, -0.1590],
            [2.00, 0.3341, 0.1260, -0.6920, 0.6877, -0.0600]
        ]
        if self.pitchType in [PitchType.TRIANGLE, PitchType.ROTATED_TRIANGLE]:
            betaTable = betaTableTriangular
        elif self.pitchType in [PitchType.SQUARE, PitchType.ROTATED_SQUARE]:
            betaTable = betaTableSquare
        else:
            raise ValueError("Invalid pitch type <<%r>> given, which is not in pitch enum PitchType" % self.pitchType)

        betaOut = None
        hpRatio = self.__hpRatio()

        for i in range(len(betaTable)):
            betaRow = betaTable[i]
            if i == 0:
                if hpRatio <= betaRow[0]:
                    betaOut = betaRow[number]

            elif i == len(betaTable) and not hpRatio < betaRow[0]:
                betaOut = betaRow

            else:
                if hpRatio < betaRow[0]:
                    betaOut = self.__interpolate(betaTable[i - 1][0], betaTable[i - 1][number], betaTable[i][0], betaTable[i][number], hpRatio)

        return betaOut

    def E_star_ratio(self):
        """

        :return: Effective modulus of elasticity for perforated region of tubesheet
        """
        alpha_0 = self.__alpha(0)
        alpha_1 = self.__alpha(1)
        alpha_2 = self.__alpha(2)
        alpha_3 = self.__alpha(3)
        alpha_4 = self.__alpha(4)
        mu_star = self.mu_star()


        return alpha_0 + alpha_1 * mu_star + alpha_2 * (mu_star ** 2) + alpha_3 * (mu_star ** 3) + alpha_4 * (mu_star ** 4)

    def E_star(self):
        return self.E_star_ratio() * self.E


    def vu_star(self):
        """

        :return: Effective poisson's ratio of perforated region of the tubesheet
        """
        beta_0 = self.__beta(0)
        beta_1 = self.__beta(1)
        beta_2 = self.__beta(2)
        beta_3 = self.__beta(3)
        beta_4 = self.__beta(4)
        mu_star = self.mu_star()

        return beta_0 + beta_1 * mu_star + beta_2 * (mu_star ** 2) + beta_3 * (mu_star ** 3) + beta_4 * (mu_star ** 4)

if __name__ == "__main__":
    params = UHX11Params(
            radius_to_outermost_tube_hole_center=30.00,
            nominal_tube_OD=1.00,
            nominal_tube_wall_thickness=0.083,
            modulus_of_elasticity_of_tubes_at_tubesheet_design_temperature=26900000,
            modulus_of_elasticity_for_tubesheet_material_at_tubesheet_design_temperature=26900000,
            allowable_stress_of_tubes_at_tubesheet_design_temperature=13400,
            allowable_stress_for_tubesheet_material_at_tubesheet_design_temperature=19000,
            tube_pitch=1.25,
            area_of_untubed_lanes=99.00,
            tube_side_pass_partition_groove_depth=0,
            tubesheet_corrosion_tube_side=0.5,
            tubesheet_thickness=5,
            expanded_depth_of_tube_in_tubesheet=(5 * 1),
            pitch_type=PitchType.SQUARE)
    calcs = UHX11(params)

    print("rho = " + str(calcs.rho()))
    print("D_o = " + str(calcs.D_o()))
    print("mu = " + str(calcs.mu()))
    print("p_star = " + str(calcs.p_star()))
    print("d_star = " + str(calcs.d_star()))
    print("mu_star = " + str(calcs.mu_star()))
    print("hg_prime = " + str(calcs.hg_prime()))
    print("h/p ratio = " + str(calcs.h / calcs.p))
    print("E_star_ratio = " + str(calcs.E_star_ratio()))
    print("E_star = " + str(calcs.E_star()))
    print("vu_star = " + str(calcs.vu_star()))