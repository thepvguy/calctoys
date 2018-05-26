import math
import enum


class SaddleAnalysisBuilder:
    def __init__(self):
        pass


class SaddleParams:
    def __init__(
            self,
            component_name,
            which_side,
            internal_pressure,
            shell_modulus_of_elasticity,
            saddle_welded_to_vessel,
            tan_dist,
            head_depth,
            head_type,
            vessel_tangent_length,
            outside_radius,
            shell_thickness,
            head_thickness,
            saddle_width,
            web_thickness,
            head_inside_radius,
            wear_plate_present,
            wear_plate_welded_to_shell,
            wear_plate_thickness,
            wear_plate_angle,
            wear_plate_width,
            base_plate_thickness,
            load,
            contact_angle_degrees,
            is_exceptional_condition,
            wear_plate_allowable,
            cylinder_allowable_stress,
            rings = None
    ):
        self.component_name = component_name
        self.which_side = which_side
        self.internal_pressure = internal_pressure
        self.shell_modulus_of_elasticity = shell_modulus_of_elasticity
        self.tan_dist = tan_dist
        self.saddle_welded_to_vessel = saddle_welded_to_vessel
        self.head_depth = head_depth
        self.head_type = head_type
        self.vessel_tangent_length = vessel_tangent_length
        self.outside_radius = outside_radius  # In inches
        self.head_thickness = head_thickness
        self.shell_thickness = shell_thickness
        self.saddle_width = saddle_width
        self.web_thickness = web_thickness
        self.head_inside_radius = head_inside_radius
        self.wear_plate_present = wear_plate_present
        self.wear_plate_welded_to_shell = wear_plate_welded_to_shell
        self.wear_plate_thickness = wear_plate_thickness
        self.wear_plate_width = wear_plate_width
        self.wear_plate_angle = wear_plate_angle
        self.base_plate_thickness = base_plate_thickness
        self.load = load  # Load on one saddle
        self.contact_angle_degrees = contact_angle_degrees
        self.is_exceptional_condition = is_exceptional_condition
        self.rings = rings
        self.wear_plate_allowable = wear_plate_allowable
        self.cylinder_allowable_stress = cylinder_allowable_stress

        if not self.isSane():
            raise ValueError("Zero or None passed into SaddleParams object")

    def isSane(self):
        isSane = True

        if self.outside_radius == 0 or self.outside_radius is None:
            isSane = False

        return isSane


class SaddleCalcsDiv2:
    def __init__(self, params: SaddleParams):
        # Design conditions
        self.P = params.internal_pressure
        self.E_y = params.shell_modulus_of_elasticity
        self.Q = params.load  # Load on one saddle
        self.isExceptional = params.is_exceptional_condition
        if params.saddle_welded_to_vessel:
            self.k = 0.1
        else:
            self.k = 1.0

        # Vessel Info
        self.which_side = params.which_side
        self.a = params.tan_dist
        self.h_2 = params.head_depth
        self.headType = params.head_type
        self.L = params.vessel_tangent_length
        self.t_h = params.head_thickness
        self.t_s = params.shell_thickness
        self.R_i = params.head_inside_radius
        self.S = params.cylinder_allowable_stress

        # Geometry info
        self.b = params.saddle_width
        self.t_w = params.web_thickness
        self.wpIsPresent = params.wear_plate_present
        self.t_basePlate = params.base_plate_thickness
        self.theta_deg = params.contact_angle_degrees
        self.rings = None  # Do ring stuff later

        if not self.isSane():
            raise ValueError("One of the inputs is zero or None")

        if self.wpIsPresent:
            self.wpIsWeldedToShell = params.wear_plate_welded_to_shell
            self.t_r = params.wear_plate_thickness
            self.b_1 = params.wear_plate_width
            self.theta_1_deg = params.wear_plate_angle
            self.theta_1 = math.radians(self.theta_1_deg)
            self.S_r = params.wear_plate_allowable
            self.beta_1 = self.calcBeta_1()
        else:
            self.wpIsWeldedToShell = None
            self.t_r = None
            self.b_1 = None
            self.theta_1_deg = None
            self.theta_1 = None
            self.beta_1 = None

        self.theta = math.radians(self.theta_deg)
        self.R_m = params.outside_radius - params.shell_thickness * 0.5
        self.alpha = self.calcAlpha()
        self.beta = self.calcBeta()
        self.delta = self.calcDelta()
        self.rho = self.calcRho()

    def isSane(self):
        isSane = True

        if self.theta_deg < 120:
            isSane = False

        return isSane

    def thetaToRadians(self):
        return math.radians(self.theta_deg)

    def calcAlpha(self):
        return 0.95 * (math.pi - 0.5 * self.theta)

    def calcBeta(self):
        return math.pi - (0.5 * self.theta)

    def calcBeta_1(self):
        return math.pi - (0.5 * self.theta_1)

    def calcDelta(self):
        return (math.pi / 6.0) + (5.0 * self.theta / 12.0)

    def calcRho(self):
        errorInterval = 0.0000001
        oldRho = 0.5 * math.pi
        rho = 0.0
        denom = 0.5 + (math.pi - self.beta) * (1 / math.tan(self.beta))
        while abs(oldRho - rho) > errorInterval:
            oldRho = rho
            temp = math.atan(rho / denom)
            if temp <= 0:
                rho = math.pi + temp
            else:
                rho = temp
        if (rho - errorInterval) > 0:
            return rho
        else:
            return errorInterval  # Zero breaks stuff

    def calcK1(self):
        num = self.delta + math.sin(self.delta) * math.cos(self.delta) - \
              ((2.0 * math.sin(self.delta) ** 2.0) / self.delta)
        den = math.pi * ((math.sin(self.delta) / self.delta) - math.cos(self.delta))
        return num / den

    def calcK1Prime(self):
        num = self.delta + math.sin(self.delta) * math.cos(self.delta) - (
        (2.0 * math.sin(self.delta) ** 2.0) / self.delta)
        den = math.pi * (1.0 - (math.sin(self.delta) / self.delta))
        return num / den

    def calcK2(self):
        return math.sin(self.alpha) / (math.pi - self.alpha + math.sin(self.alpha) * math.cos(self.alpha))

    def calcK3(self):
        return (math.sin(self.alpha) / math.pi) * \
               ((self.alpha - math.sin(self.alpha) * math.cos(self.alpha)) /
                (math.pi - self.alpha + math.sin(self.alpha) * math.cos(self.alpha)))

    def calcK4(self):
        return 0.375 * (math.sin(self.alpha) ** 2) / (
        math.pi - self.alpha + math.sin(self.alpha) * math.cos(self.alpha))

    def calcK5(self):
        return (1.0 + math.cos(self.alpha)) / (math.pi - self.alpha + math.sin(self.alpha) * math.cos(self.alpha))

    def calcK6(self):

        sinBeta = math.sin(self.beta)
        cosBeta = math.cos(self.beta)
        sinBetaOverBeta = sinBeta / self.beta
        sin2BetaOverBeta = math.sin(2.0 * self.beta) / self.beta
        cosBetaOverBeta = cosBeta / self.beta

        term1 = 0.75 * cosBeta * sinBetaOverBeta ** 2.0
        term2 = 1.25 * sinBeta * cosBeta * cosBetaOverBeta
        term3 = 0.5 * cosBeta ** 3.0
        term4 = 0.25 * sinBetaOverBeta
        term5 = 0.25 * cosBeta
        term6 = self.beta * sinBeta
        term7 = (sinBetaOverBeta ** 2.0) - 0.5 - 0.25 * sin2BetaOverBeta
        num = term1 - term2 + term3 - term4 + term5 - term6 * term7
        den = 2.0 * math.pi * term7
        return num / den

    def calcK6_1(self):

        sinBeta = math.sin(self.beta_1)
        cosBeta = math.cos(self.beta_1)
        sinBetaOverBeta = sinBeta / self.beta_1
        sin2BetaOverBeta = math.sin(2.0 * self.beta_1) / self.beta_1
        cosBetaOverBeta = cosBeta / self.beta_1

        term1 = 0.75 * cosBeta * sinBetaOverBeta ** 2.0
        term2 = 1.25 * sinBeta * cosBeta * cosBetaOverBeta
        term3 = 0.5 * cosBeta ** 3.0
        term4 = 0.25 * sinBetaOverBeta
        term5 = 0.25 * cosBeta
        term6 = self.beta_1 * sinBeta
        term7 = (sinBetaOverBeta ** 2.0) - 0.5 - 0.25 * sin2BetaOverBeta
        num = term1 - term2 + term3 - term4 + term5 - term6 * term7
        den = 2.0 * math.pi * term7
        return num / den

    def calcK7(self):

        tanRadRatio = self.a / self.R_m

        if tanRadRatio <= 0.5:
            return self.calcK6() * 0.25

        elif tanRadRatio >= 1:
            return self.calcK6()

        else:  # 0.5 < tanRadRatio < 1
            return 1.5 * self.calcK6() * tanRadRatio - 0.5 * self.calcK6()

    def calcK7_1(self):
        """

        :return: K7 evaluated at the wear plate angle
        """
        tanRadRatio = self.a / self.R_m

        if tanRadRatio <= 0.5:
            return self.calcK6_1() * 0.25

        elif tanRadRatio >= 1.0:
            return self.calcK6_1()

        else:  # 0.5 < tanRadRatio < 1
            return 1.5 * self.calcK6_1() * tanRadRatio - 0.5 * self.calcK6_1()

    def calcK8(self):
        num = math.cos(self.beta) * \
              (1.0 -
               0.25 * math.cos(2.0 * self.beta) +
               (2.25 * math.sin(self.beta) * math.cos(self.beta) / self.beta) -
               3.0 * ((math.sin(self.beta) / self.beta) ** 2))
        den = 2.0 * math.pi * (((math.sin(self.beta) / self.beta) ** 2) -
                               0.5 - (math.sin(2.0 * self.beta) / 4.0 * self.beta))

        return (num / den) + (self.beta * math.sin(self.beta) / (2.0 * math.pi))

    def calcK9(self):  # woof woof
        return (1.0 / (2.0 * math.pi)) * \
               ((-0.5 + ((math.pi - self.beta) / math.tan(self.beta))) *
                math.cos(self.rho) - (self.rho * math.sin(self.rho)))

    def calcK10(self):
        return (1.0 / (2.0 * math.pi)) * (self.rho * math.sin(self.rho) + math.cos(self.rho) *
                                      (1.5 + (math.pi - self.beta) * (1.0 / math.tan(self.beta))) -
                                      ((math.pi - self.beta) / (math.sin(self.beta))))

    def calc_b1Min(self):
        """

        :return: minimum b_1 per equation 4.15.1
        """
        return min(self.b + 1.56 * math.sqrt(self.R_m * self.t_s), 2.0 * self.a)

    def calcTheta1Min(self):
        """

        :return: minimum included angle of the wear plate, if present per equation 4.15.2
        """
        return self.theta + (self.theta / 12.0)

    def calcM1(self):
        """

        :return: Moment at the saddle per equation 4.15.3
        """
        return -self.Q * self.a * (
        1.0 - (1.0 - (self.a / self.L) + (((self.R_m ** 2) - (self.h_2 ** 2)) / (2.0 * self.a * self.L))) / (
        1.0 + ((4.0 * self.h_2) / (3.0 * self.L))))

    def calcM2(self):
        """

        :return: Moment in the middle of the vessel per equation 4.15.4
        """
        return 0.25 * self.Q * self.L * (
        ((1.0 + ((2.0 * ((self.R_m ** 2) - (self.h_2 ** 2))) / (self.L ** 2))) / (1.0 + ((4.0 * self.h_2) / (3.0 * self.L)))) - (
        4.0 * self.a / self.L))

    def calcT(self):
        """

        :return: Shear force at the saddle per equation 4.15.5
        """
        return (self.Q * (self.L - 2.0 * self.a)) / (self.L + ((4.0 * self.h_2) / 3.0))

    def calcSigma1(self):
        """

        :return: Longitudinal stress at the top of the shell between the saddles per equation 4.15.6
        """
        return (self.P * self.R_m / (2.0 * self.t_s)) - (self.calcM2() / (math.pi * (self.R_m ** 2) * self.t_s))

    def calcSigma2(self):
        """

        :return: Longitudinal Stress at the bottom of the shell between the saddles per equation 4.15.7
        """

        return (self.P * self.R_m / (2.0 * self.t_s)) + (self.calcM2() / (math.pi * (self.R_m ** 2) * self.t_s))

    def calcSigma3(self):
        """

        :return: Longitudinal stress at the top of a stiffened shell at the saddle per equation 4.15.8
        """

        return (self.P * self.R_m / (2.0 * self.t_s)) - (self.calcM1() / (math.pi * (self.R_m ** 2) * self.t_s))

    def calcSigma4(self):
        """

        :return: Longitudinal stress at the bottom of a stiffened shell at the saddle per equation 4.15.9
        """

        return (self.P * self.R_m / (2.0 * self.t_s)) + (self.calcM1() / (math.pi * (self.R_m ** 2) * self.t_s))

    def calcSigma3Star(self):
        """

        :return: Longitudinal stress at the top of an unstiffened shell at the saddle per equation 4.15.10
        """
        return (self.P * self.R_m / (2.0 * self.t_s)) - \
               (self.calcM1() / (self.calcK1() * math.pi * (self.R_m ** 2) * self.t_s))

    def calcSigma4Star(self):
        """

        :return: Longitudinal stress at the bottom of an unstiffened shell at the saddle per equation 4.15.11
        """
        return (self.P * self.R_m / (2.0 * self.t_s)) + (
        self.calcM1() / (self.calcK1Prime() * math.pi * (self.R_m ** 2) * self.t_s))

    def calcSc(self):
        """

        :return: Max compressive stress  in the shell per equation 4.15.12
        """
        if self.isExceptional:
            K = 1.35
        else:
            K = 1.0

        return (K * self.t_s * self.E_y) / (16.0 * self.R_m)

    def calcTau1(self):
        """

        :return: Max shear stress in a shell with a stiffening ring in the plane of the saddle per equation 4.15.13
        """
        return self.calcT() / (math.pi * self.R_m * self.t_s)

    def calcTau2(self):
        """

        :return: Max shear stress in a shell stiffened by two rings, or not stiffened at all, per equation 4.15.14
        """
        return (self.calcK2() * self.calcT()) / (self.R_m * self.t_s)

    def calcTau3(self):
        """

        :return: Max shear stress in the shell stiffened by a head or tubesheet per equation 4.15.16
        """
        return (self.calcK3() * self.Q) / (self.R_m * self.t_s)

    def calcTau3Star(self):
        """

        :return: Extra shear stress in a head or tubesheetstiffening a shell per equation 4.15.17
        """
        return (self.calcK3() * self.Q) / (self.R_m * self.t_h)

    def calcSigma5(self):
        """

        :return: Membrane stress in a torispherical, elliptical, or flat head per equation 4.15.17, 18, or 19
        """
        if self.headType is HeadType.TORISPHERICAL:
            return ((self.calcK4() * self.Q) / (self.R_m * self.t_h)) + ((self.P * self.R_i) / (2.0 * self.t_h))

        elif self.headType in [HeadType.ELLIPTICAL, HeadType.PIPE_CAP, HeadType.HEMISPHERICAL]:
            return ((self.calcK4() * self.Q) / (self.R_m * self.t_h)) + (
            (self.P * self.R_i ** 2) / (2.0 * self.t_h * self.h_2))

        elif self.headType is HeadType.FLAT_COVER or self.headType is HeadType.TUBESHEET:
            return 0.0

        else:
            raise ValueError("Value %r is not an acceptable head type." % self.headType)

    def calcMBeta(self):
        """

        :return: MAx circumferential bending moment at the saddle per equation 4.15.20 or 21
        """
        if self.rings is None or self.rings.count == 0 or self.rings.count == 1:
            return self.calcK7() * self.Q * self.R_m
        elif self.rings.count > 1.0:
            return self.calcK10() * self.Q * self.R_m
        else:
            raise ValueError("Cannot have %r quantity of rings" % self.rings.count)

    def calcX1(self):
        # TODO: fine tune this
        """

        :return: Contributing length of the shell adjacent to the left saddle per equation 4.15.22
        """
        return 0.78 * math.sqrt(self.R_m * self.t_s)

    def calcX2(self):
        # TODO: fine tune this
        """

        :return: Contributing length of the shell adjacent to the right saddle per equation 4.15.22
        """
        return 0.78 * math.sqrt(self.R_m * self.t_s)

    def calcSigma6(self):
        """

        :return: max compressive circumferential membrane stress
        in the cylindrical shell at the base of the saddle per equation 4.15.23,
        or equation 4.15.37 if there are two rings.
        """

        if self.rings is not None and self.rings.count == 2:
            return (-self.calcK5() * self.Q * self.k) / (self.t_s * (self.b + 2.0 * self.calcX2()))
        else:
            return (-self.calcK5() * self.Q * self.k) / (self.t_s * (self.b + self.calcX1() + self.calcX2()))

    def calcSigma7(self):
        """

        :return: Circumferential compressive membrane stress plus bending at saddle horns
        in the plane of the saddle per equation 4.15.24 (when the vessel length is greater than 8 R_m)
        """
        return -(self.Q / (4.0 * self.t_s * (self.b + self.calcX1() + self.calcX2()))) - \
               ((3.0 * self.calcK7() * self.Q) / (2.0 * self.t_s ** 2))

    def calcSigma7Star(self):
        """

        :return: Circumferential compressive membrane stress plus bending at saddle horns
        in the plane of the saddle per equation 4.15.25 (when the vessel length is less than 8 R_m)
        """
        return -(self.Q / (4.0 * self.t_s * (self.b + self.calcX1() + self.calcX2()))) - \
               ((12.0 * self.calcK7() * self.Q * self.R_m) / (self.L * self.t_s ** 2))

    def calcEta(self):
        """

        :return: Ratio between the reinforcing plate allowable stress at temperature
        and the cylindrical shell material, or one, whichever is less, per equation 4.15.29
        """
        return min((self.S_r / self.S), 1.0)

    def calcSigma6r(self):
        """

        :return: Max compressive circumferential membrane stress in the cylindrical shell
        with a wear plate at the base of the saddle per equation 4.15.26
        """
        return (- self.calcK5() * self.Q * self.k) / (self.b_1 * (self.t_s + self.calcEta() * self.t_r))

    def calcSigma7r(self):
        """

        :return: Circumferential compressive membrane stress plus bending at saddle horns
        in the plane of the saddle when the vessel has a wear plate
         per equation 4.15.27  (when the vessel length is greater than 8 R_m)
        """
        return (-self.Q / (4.0 * (self.t_s + self.calcEta() * self.t_r) * self.b_1)) - \
               ((3.0 * self.calcK7() * self.Q) / (2.0 * (self.t_s + self.calcEta() * self.t_r) ** 2))

    def calcSigma7rStar(self):
        """

        :return:Circumferential compressive membrane stress plus bending at saddle horns
        in the plane of the saddle when the vessel has a wear plate
        per equation 4.15.28 (when the vessel length is less than 8 R_m)
        """
        return (-self.Q / (4.0 * (self.t_s + self.calcEta() * self.t_r) * self.b_1)) - \
               ((12.0 * self.calcK7() * self.Q * self.R_m) / (self.L * (self.t_s + self.calcEta() * self.t_r) ** 2))

    def calcSigma7_1(self):
        """

        :return: Circumferential compressive membrane plus bending stress
        at the ends of the reinforcing plates in a vessel without stiffening rings
         where the vessel is less than 4 diameters long, per equation 4.15.30
        """
        return ((- self.Q) / (4.0 * self.t_s * (self.b + self.calcX1() + self.calcX2()))) - \
               ((3.0 * self.calcK7_1() * self.Q) / (2.0 * self.t_s ** 2))

    def calcSigma7Star_1(self):
        """

        :return: Circumferential compressive membrane plus bending stress
        at the ends of the reinforcing plates in a vessel without stiffening rings
        where the vessel is greater than 4 diameters long, per equation 4.15.31
        """
        return ((- self.Q) / (4.0 * self.t_s * (self.b + self.calcX1() + self.calcX2()))) - \
               ((12.0 * self.calcK7_1() * self.Q * self.R_m) / (self.L * self.t_s ** 2))

    def calcSigma6Star(self):
        """

        :return: Circumferential compressive membrane plus bending stress
        along the plane of the saddle support in the shell with a ring in the plane of the saddle per equation 4.15.32
        """
        return (- self.calcK5() * self.Q * self.k) / self.rings.A

    def calcSigma8(self):
        """

        :return: The Circumferential compressive membrane plus bending stress in the shell
        for stiffening rings inside the vessel per equation 4.15.33
        """
        return (-self.calcK8() * self.Q / self.rings.A) - \
               ((self.calcK6() * self.Q * self.R_m * self.rings.c_1) / self.rings.I)

    def calcSigma9(self):
        """

        :return: The Circumferential compressive membrane plus bending stress in the rings
        for stiffening rings inside the vessel per equation 4.15.34
        """
        return (-self.calcK8() * self.Q / self.rings.A) + \
               ((self.calcK6() * self.Q * self.R_m * self.rings.c_2) / self.rings.I)

    def calcSigma8Star(self):
        """

        :return: The Circumferential compressive membrane plus bending stress in the shell
        for stiffening rings outside the vessel per equation 4.15.35
        """
        return (-self.calcK8() * self.Q / self.rings.A) + \
               ((self.calcK6() * self.Q * self.R_m * self.rings.c_1) / self.rings.I)

    def calcSigma9Star(self):
        """

        :return: The Circumferential compressive membrane plus bending stress in the rings
        for stiffening rings outside the vessel per equation 4.15.36
        """
        return (-self.calcK8() * self.Q / self.rings.A) - \
               ((self.calcK6() * self.Q * self.R_m * self.rings.c_2) / self.rings.I)

    def calcSigma10(self):
        """

        :return: Circumferential compressive membrane plus bending stress in the shell at xxx for stiffening rings
        located on the inside of the shell per equation 4.15.38
        """
        return (-self.calcK9() * self.Q / self.rings.A) + \
               ((self.calcK10() * self.Q * self.R_m * self.rings.c_1) / self.rings.I)

    def calcSigma11(self):
        """

        :return: Circumferential compressive membrane plus bending stress in the ring at xxx for stiffening rings
        located on the inside of the shell per equation 4.15.39
        """
        return (- self.calcK9() * self.Q / self.rings.A) - \
               ((self.calcK10() * self.Q * self.R_m * self.rings.c_2) / self.rings.I)

    def calcSigma10Star(self):
        """

        :return: Circumferential compressive membrane plus bending stress in the shell at xxx for stiffening rings
        located on the inside of the shell per equation 4.15.40
        """
        return (-self.calcK9() * self.Q / self.rings.A) - \
               ((self.calcK10() * self.Q * self.R_m * self.rings.c_1) / self.rings.I)

    def calcSigma11Star(self):
        """

        :return: Circumferential compressive membrane plus bending stress in the shell at xxx for stiffening rings
        located on the inside of the shell per equation 4.15.41
        """
        return (- self.calcK9() * self.Q / self.rings.A) + \
               ((self.calcK10() * self.Q * self.R_m * self.rings.c_2) / self.rings.I)

    def calcF_h(self):
        """

        :return: The horizontal force at the minimum section at the low point of the saddle per 4.15.42
        """
        return self.Q * ((1.0 + math.cos(self.beta) - 0.5 * math.sin(self.beta) ** 2) /
                         (math.pi - self.beta + math.sin(self.beta) * math.cos(self.beta)))

    def buildReportContext(self):
        """

        :return: Returns a structure that can be rendered to a report. TODO: more documentation here.
        """

        #TODO: Clean up the control flow. It follows the same structure as the ASME code currently.
        twoRingsCountedAsOne = False
        if self.rings is not None and self.rings.count == 2 and self.rings.h < 1.56 * math.sqrt(self.R_m * self.t_s):
            twoRingsCountedAsOne = True

        result = {}

        forces = {
            "M1": self.calcM1(),
            "M2": self.calcM2(),
            "T": self.calcT(),
            "F_h": self.calcF_h()
        }

        stresses = {
            "Sigma1": self.calcSigma1(),
            "Sigma2": self.calcSigma2()
        }

        constants = {
            "Theta": self.thetaToRadians(),
            "Alpha": self.calcAlpha(),
            "Beta": self.calcBeta(),
            "Delta": self.calcDelta(),
            "Rho": self.calcRho()
        }

        if self.rings is not None or self.a <= 0.5 * self.R_m:
            stresses["Sigma3"] = self.calcSigma3()
            stresses["Sigma4"] = self.calcSigma4()
        else:
            constants["K1"] = self.calcK1() # TODO: Make sure this is right, Pretty sure it's supposed to go with S1 and S2
            constants["K1Prime"] = self.calcK1Prime()
            stresses["Sigma3Star"] = self.calcSigma3Star()
            stresses["Sigma4Star"] = self.calcSigma4Star()

        if self.rings is not None and self.rings.count != 0:
            if self.rings.count == 1 or twoRingsCountedAsOne:
                stresses["Tau1"] = self.calcTau1()
            elif self.rings.count == 2:
                constants["K2"] = self.calcK2()
                stresses["Tau2"] = self.calcTau2()
            else:
                raise ValueError("Invalid ring configuration")

        else:
            if self.a > 0.5*self.R_m:
                stresses["Tau2"] = self.calcTau2()
            else:
                constants["K3"] = self.calcK3()
                constants["K4"] = self.calcK4()
                stresses["Tau3"] = self.calcTau3()
                stresses["Tau3Star"] = self.calcTau3Star()
                stresses["Sigma5"] = self.calcSigma5()

        if self.rings is None or self.rings.count == 1 or twoRingsCountedAsOne:
            constants["K7"] = self.calcK7()
        else:
            constants["K10"] = self.calcK10()

        stresses["MBeta"] = self.calcMBeta()

        if self.rings is None or self.rings.count == 0:

            if self.wpIsPresent and self.wpIsWeldedToShell:
                constants["Eta"] = self.calcEta()
                constants["K6"] = self.calcK6()
                stresses["Sigma6r"] = self.calcSigma6r()
                if self.t_r > self.t_s * 2:
                    constants["K7_1"] = self.calcK7_1()
                    if self.L > 8 * self.R_m:
                        stresses["Sigma7_1"] = self.calcSigma7_1()
                    else:
                        stresses["Sigma7Star_1"] = self.calcSigma7Star_1()
                else:
                    if self.L >= 8 * self.R_m:
                        stresses["Sigma7r"] = self.calcSigma7r()
                    else:
                        stresses["Sigma7rStar"] = self.calcSigma7rStar()
            else:
                constants["K6"] = self.calcK6()
                stresses["Sigma6"] = self.calcSigma6()

                if self.L >= 8 * self.R_m:
                    stresses["Sigma7"] = self.calcSigma7()
                else:
                    stresses["Sigma7Star"] = self.calcSigma7Star()

        elif self.rings.count == 1 or twoRingsCountedAsOne:
            constants["K6"] = self.calcK6()
            stresses["Sigma6Star"] = self.calcSigma6Star()
            constants["K8"] = self.calcK8()
            if self.rings.isInternal:
                stresses["Sigma8"] = self.calcSigma8()
                stresses["Sigma9"] = self.calcSigma9()

            elif not self.rings.isInternal:
                stresses["Sigma8Star"] = self.calcSigma8Star()
                stresses["Sigma9Star"] = self.calcSigma9Star()

            else:
                raise ValueError("Rings should be either internal or external")

        elif self.rings.count == 2:
            constants["K6"] = self.calcK6()
            stresses["Sigma6"] = self.calcSigma6()
            constants["K9"] = self.calcK9()
            if self.rings.isInternal:
                stresses["Sigma10"] = self.calcSigma10()
                stresses["Sigma11"] = self.calcSigma11()

            elif not self.rings.isInternal:
                stresses["Sigma10Star"] = self.calcSigma10Star()
                stresses["Sigma11Star"] = self.calcSigma11Star()

            else:
                raise ValueError("Rings should be either internal or external")

        else:
            raise ValueError("Invalid ring configuration location 2")

        acceptance = {
            "_b1Min": self.calc_b1Min(),
            "Theta1Min": self.calcTheta1Min(),
            "Sc": self.calcSc(),
            "X1": self.calcX1(),
            "X2": self.calcX2()
        }

        result['constants'] = constants
        result['forces'] = forces
        result['stresses'] = stresses
        result['acceptance'] = acceptance

        return result


class SaddleRings:
    """
    Class for saddle ring data.
    """

    def __init__(self,
                 count,
                 cross_sectional_area_one_ring,
                 dist_from_shell_to_ring_neutral_axis,
                 dist_from_neutral_axis_to_ring_outer_edge,
                 area_moment_of_inertia,
                 allowable_stress,
                 section_modulus,
                 spacing_center_to_center=0,
                 is_internal=False
                 ):
        self.count = count
        self.S_s = allowable_stress
        self.A = cross_sectional_area_one_ring
        self.c_1 = dist_from_shell_to_ring_neutral_axis
        self.c_2 = dist_from_neutral_axis_to_ring_outer_edge
        self.I = area_moment_of_inertia
        self.Z = section_modulus  # one ring
        self.isInternal = is_internal

        if self.count > 1:
            self.h = spacing_center_to_center
        else:
            self.h = None


class HeadType(enum.Enum):
    ELLIPTICAL = 1
    HEMISPHERICAL = 2
    TORISPHERICAL = 3
    FLAT_COVER = 4
    PIPE_CAP = 5
    TUBESHEET = 100


class Side(enum.Enum):
    RIGHT = 1
    LEFT = 2


class SaddleReportStrings(enum.Enum):
    A = "A"
    a = "a"
    b = "b"
    b_1 = "b_1"
    c_1 = "c_1"
    c_2 = "c_2"
    E_y = "E_y"
    E = "E"
    eta = "eta"
    F_h = "F_h"
    h = "h"
    h_2 = "h_2"
    I = "I"
    k = "k"
    K = "K"
    L = "L"
    M_1 = "M_1"
    M_2 = "M_2"
    P = "P"
    Q = "Q"
    R_i = "R_i"
    R_m = "R_m"
    S = "S"
    S_c = "S_c"
    S_h = "S_h"
    S_r = "S_r"
    S_s = "S_s"
    t = "t"
    t_h = "t_h"
    t_r = "t_r"
    T = "T"
    theta = "theta"
    theta_1 = "theta_1"
    x_1 = "x_1"
    x_2 = "x_2"

    def __repr__(self):
        return '<%s, %s>' % (self.__class__.__name__, self.name)

if __name__ == "__main__":

    # Should eventually just be saddle info, vessel info will be abstracted out.
    inputs = SaddleParams(
        component_name="Saddle",
        which_side=Side.RIGHT,
        internal_pressure=386.3,
        shell_modulus_of_elasticity=24920000,
        saddle_welded_to_vessel=True,
        tan_dist=12,
        head_depth=(2/3) * 6.25,
        head_type=HeadType.ELLIPTICAL,
        vessel_tangent_length=212.3333,
        outside_radius=24.5,
        shell_thickness=0.5,
        head_thickness=0.25,
        saddle_width=10,
        web_thickness=0.5,
        head_inside_radius=12,
        wear_plate_present=True,
        wear_plate_welded_to_shell=True,
        wear_plate_thickness=0.5,
        wear_plate_angle=165,
        wear_plate_width=12,
        base_plate_thickness=0.5,
        load=1026,
        contact_angle_degrees=132,
        is_exceptional_condition=False,
        wear_plate_allowable=18800,
        cylinder_allowable_stress=18800,
        rings=None
    )

    #Should eventually be as follows:
    # SaddleCalcsDiv2(vessel, saddle_params)
    calcs = SaddleCalcsDiv2(inputs)

    print(calcs.buildReportContext())
