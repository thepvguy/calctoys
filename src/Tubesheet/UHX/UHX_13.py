from Table_13_1_and_2 import Table_13_1, Table_13_2
from _UHX_common import Configuration as cfg
import math


class UHX_13_5_calcs:
    def __init__(self, params):
        """
            Determine D_o, mu, mu_star, h_prime_g from UHX-11.5.1;
            tubesheet parameters
        """
        self.D_o = params.D_o
        self.mu = params.mu
        self.mu_star = params.mu_star
        self.__h_prime_g = params.h_prime_g
        self.vu_star = params.vu_star
        self.E = params.E
        self.E_star = params.E_star
        self.h = params.h
        self.A = params.A
        self.C = params.C
        self.S = params.S
        self.S_PS = params.S_PS
        self.W_star = params.W_star
        """
            General parameters
        """
        self.config = params.config
        self.is_operating_case = params.is_operating_case

        """
            Design Conditions
        """
        self.P_s = params.P_s
        self.P_t = params.P_t
        self.T_t_m = params.T_t_m
        self.T_a = params.T_a
        self.T_s_m = params.T_s_m

        """
            Mating Flange Parameters; can possibly be None
        """
        self.G_c = params.G_c
        self.G_s = params.G_s
        self.G_1 = params.G_1

        """
            Adjacent shell parameters
        """
        self.D_s = params.D_s
        self.t_s = params.t_s
        self.E_s = params.E_s
        self.L = params.L
        self.vu_s = params.vu_s
        self.alpha_s_m = params.alpha_s_m
        self.S_PS_s = params.S_PS_s
        self.E_s_w = params.E_s_w
        self.S_s_b = params.S_s_b
        self.S_s = params.S_s

        """
            Adjacent channel parameters
        """
        self.D_c = params.D_c
        self.vu_c = params.vu_c
        self.t_c = params.t_c
        self.E_c = params.E_c
        self.S_c = params.S_c
        self.S_PS_c = params.S_PS_c

        """
            Tube parameters
        """
        self.N_t = params.N_t
        self.d_t = params.d_t
        self.t_t = params.t_t
        self.E_t = params.E_t
        self.alpha_t_m = params.alpha_t_m
        self.vu_t = params.vu_t
        self.A_p = params.A_p
        self.C_p = params.C_p
        self.L_max = params.L_max
        self.S_t = params.S_t
        self.k = params.k
        self.l = params.l
        self.S_y_t = params.S_y_t
        """
            Expansion Joint Properties
        """
        self.has_expansion_joint = params.has_expansion_joint
        self.K_J = params.K_J
        if self.has_expansion_joint is False:
            self.D_J = self.D_s
        else:
            self.D_J = params.D_J

        self.sanity_check()

        """
            Private properties
        """
        self.__Z_d = None
        self.__Z_v = None
        self.__Z_w = None
        self.__Z_m = None
        self.__F_m = None
        self.__F_t_min = None
        self.__F_t_max = None


    def sanity_check(self):
        if self.config not in [cfg.a, cfg.b, cfg.c, cfg.d]:
            raise ValueError("Improper value '%s' for parameter 'config'" % self.config)

    @property
    def h_prime_g(self):
        """
            UHX 13.5.1
        """
        if self.is_operating_case:
            return 0
        else:
            return self.__h_prime_g

    @property
    def a_o(self):
        """
        UHX 13.5.1
        :return: Equivalent Radius of the outer tube limit circle
        """
        return self.D_o / 2.0

    @property
    def a_s(self):
        """
        UHX 13.5.1
        :return: Radial shell dimension
        """
        if self.config == cfg.d:
            return 0.5 * self.G_s
        else:
            return self.D_s * 0.5

    @property
    def a_c(self):
        """
        UHX 13.5.1
        :return: Radial shell dimension
        """
        if self.config != cfg.a:
            return 0.5 * self.G_c
        else:
            return self.D_c * 0.5


    @property
    def rho_s(self):
        """
        UHX 13.5.1
        :return: Shell ID to OTL ratio
        """
        return self.a_s / self.a_o

    @property
    def rho_c(self):
        """
        UHX 13.5.1
        :return: Channel ID to OTL ratio
        """
        return self.a_c / self.a_o

    @property
    def x_s(self):
        """
        UHX 13.5.1
        :return: Ratio of shell side tube area to area inside the OTL
        """
        return 1 - self.N_t * (self.d_t / (2 * self.a_o)) ** 2

    @property
    def x_t(self):
        """
        UHX 13.5.1
        :return: Ratio of the tube side area to the area inside the OTL
        """
        return 1 - self.N_t * ((self.d_t - 2 * self.t_t)/(2*self.a_o)) ** 2

    def __step_1_results(self):
        # UHX 13.5.1
        return {
                "a_o": self.a_o,
                "a_c": self.a_c,
                "a_s": self.a_s,
                "rho_s": self.rho_s,
                "rho_c": self.rho_c,
                "x_s": self.x_s,
                "x_t": self.x_t
            }

    @property
    def K_s(self):
        """
        UHX 13.5.2
        :return: Shell stiffness
        """
        return (math.pi * self.t_s * (self.D_s + self.t_s) * self.E_s) / self.L

    @property
    def K_t(self):
        """
        UHX 13.5.2
        :return: Single tube stiffness
        """
        return (math.pi * self.t_t * (self.d_t + self.t_t) * self.E_t) / self.L

    @property
    def K_s_t(self):
        """
        UHX 13.5.2
        :return: Ratio of the shell to the tube bundle
        """
        return self.K_s / (self.N_t * self.K_t)

    @property
    def J(self):
        """
        UHX 13.5.2
        :return: effective spring rate of the expansion joint plus the shell
        """
        if not self.has_expansion_joint:
            return 1
        else:
            return 1 / (1 + (self.K_s / self.K_J))

    @property
    def beta_s(self):
        """
        UHX 13.5.2
        :return: Shell coefficient
        """
        if self.config == cfg.d:
            return 0
        else:
            return ((12 * (1 - (self.vu_s ** 2))) ** 0.25) / (((self.D_s + self.t_s)*self.t_s) ** 0.5)

    @property
    def k_s(self):
        """
        UHX 13.5.2
        :return: Shell coefficient
        """
        if self.config == cfg.d:
            return 0
        else:
            return self.beta_s * ((self.E_s * self.t_s ** 3) / (6 * (1 - self.vu_s ** 2)))

    @property
    def lambda_s(self):
        """
        UHX 13.5.2
        :return: Shell coefficient
        """
        if self.config == cfg.d:
            return 0
        else:
            return (6 * self.D_s / self.h ** 3) * self.k_s * (1 + self.h * self.beta_s + 0.5 * (self.h * self.beta_s) ** 2)

    @property
    def delta_s(self):
        """
        UHX 13.5.2
        :return: Shell coefficient
        """
        if self.config == cfg.d:
            return 0
        else:
            return ((self.D_s ** 2)/(4 * self.E_s * self.t_s)) * (1 - 0.5 * self.vu_s)


    @property
    def beta_c(self):
        """
        UHX 13.5.2
        :return: Channel coefficient
        """
        if self.config != cfg.a:
            return 0
        else:
            return ((12 * (1 - (self.vu_c ** 2))) ** 0.25) / (((self.D_c + self.t_c)*self.t_c) ** 0.5)

    @property
    def k_c(self):
        """
        UHX 13.5.2
        :return: Channel coefficient
        """
        if self.config != cfg.a:
            return 0
        else:
            return self.beta_c * ((self.E_c * self.t_c ** 3) / (6 * (1 - self.vu_c ** 2)))

    @property
    def lambda_c(self):
        """
        UHX 13.5.2
        :return: Channel coefficient
        """
        if self.config != cfg.a:
            return 0
        else:
            return (6 * self.D_c / self.h ** 3) * self.k_c * (1 + self.h * self.beta_c + 0.5 * (self.h * self.beta_c) ** 2)

    @property
    def delta_c(self):
        """
        UHX 13.5.2
        :return: Channel coefficient
        """
        # TODO implement hemi head
        if self.config != cfg.a:
            return 0
        else:
            return ((self.D_c ** 2)/(4 * self.E_c * self.t_c)) * (1 - 0.5 * self.vu_c)

    def __step_2_results(self):
        return {
            "K_s": self.K_s,
            "K_t": self.K_t,
            "K_s_t": self.K_s_t,
            "J": self.J,
            "beta_s": self.beta_s,
            "k_s": self.k_s,
            "lambda_s": self.lambda_s,
            "delta_s": self.delta_s,
            "beta_c": self.beta_c,
            "k_c": self.k_c,
            "lambda_c": self.lambda_c,
            "delta_c": self.delta_c
        }

    @property
    def X_a(self):
        """
        UHX 13.5.3
        :return: X_a
        """
        return (24 * (1 - self.vu_star) * self.N_t * ((self.E_t* self.t_t * (self.d_t - self.t_t) * self.a_o ** 2)/(self.E_star * self.L * self.h ** 3))) ** 0.25

    def __gettable131values_step3(self):
        table = Table_13_1(self.X_a, self.vu_star)

        self.__Z_d = float(table.Z_d)
        self.__Z_m = float(table.Z_m)
        self.__Z_v = float(table.Z_v)
        self.__Z_w = float(table.Z_w)


    @property
    def Z_d(self):
        """
        UHX 13.5.3
        :return: Z_d
        """
        if self.__Z_d is None:
            self.__gettable131values_step3()

        return self.__Z_d

    @property
    def Z_m(self):
        """
        UHX 13.5.3
        :return: Z_m
        """
        if self.__Z_m is None:
            self.__gettable131values_step3()

        return self.__Z_m

    @property
    def Z_v(self):
        """
        UHX 13.5.3
        :return: Z_v
        """
        if self.__Z_v is None:
            self.__gettable131values_step3()

        return self.__Z_v

    @property
    def Z_w(self):
        """
        UHX 13.5.3
        :return: Z_w
        """
        if self.__Z_w is None:
            self.__gettable131values_step3()

        return self.__Z_w

    def __step_3_results(self):
        return {
            "X_a": self.X_a,
            "Z_d": self.Z_d,
            "Z_m": self.Z_m,
            "Z_v": self.Z_v,
            "Z_w": self.Z_w
        }

    @property
    def K(self):
        """
        UHX 13.5.4
        :return: Diameter Ratio
        """
        return self.A / self.D_o

    @property
    def F(self):
        """
        UHX 13.5.4
        :return: Coefficient F
        """
        return ((1 - self.vu_star) / self.E_star) * (self.lambda_s + self.lambda_c + self.E * math.log(self.K, math.e))

    @property
    def phi(self):
        """
        UHX 13.5.4
        :return: Phi
        """
        return (1 - self.vu_star) * self.F

    @property
    def Q_1(self):
        """
        UHX 13.5.4
        :return: Q_1
        """
        return (self.rho_s - 1 - self.phi * self.Z_v) / (1 + self.phi * self.Z_m)

    @property
    def Q_Z_1(self):
        """
        UHX 13.5.4
        :return: Q_Z_!
        """
        return 0.5 * (self.Z_d + self.Q_1 * self.Z_w) * self.X_a ** 4

    @property
    def Q_Z_2(self):
        """
        UHX 13.5.4
        :return: Q_Z_2
        """
        return 0.5 * (self.Z_v + self.Q_1 * self.Z_m) * self.X_a ** 4

    @property
    def U(self):
        """
        UHX 13.5.4
        :return: U
        """
        return ((self.Z_w + (self.rho_s - 1) * self.Z_m) * self.X_a ** 4) / (1 + self.phi * self.Z_m)

    def __step_4_results(self):
        return {
            "K": self.K,
            "F": self.F,
            "phi": self.phi,
            "Q_1": self.Q_1,
            "Q_Z_1": self.Q_Z_1,
            "U": self.U
        }

    @property
    def gamma(self):
        """
        UHX 13.5.5
        :return: Differential axial thermal expansion
        """
        if not self.is_operating_case:
            return 0
        else:
            return (self.alpha_t_m * (self.T_t_m - self.T_a) - self.alpha_s_m*(self.T_s_m - self.T_a)) * self.L

    @property
    def omega_s(self):
        """
        UHX 13.5.5
        :return: omega_s
        """
        return self.rho_s * self.k_s * self.beta_s * self.delta_s * (1 - self.h * self.beta_s)

    @property
    def omega_s_star(self):
        """
        UHX 13.5.5
        :return: omega_s_star
        """
        return 0.25 * (self.a_o ** 2) * ((self.rho_s ** 2) - 1) * (self.rho_s - 1) - self.omega_s

    @property
    def omega_c(self):
        """
        UHX 13.5.5
        :return: omega_c
        """
        return self.rho_c * self.k_c * self.beta_c * (1 + self.beta_c * self.h)

    @property
    def omega_c_star(self):
        """
        UHX 13.5.5
        :return: omega_c_star
        """
        return (self.a_o ** 2) * ((0.25 * ((self.rho_c ** 2) + 1) * (self.rho_c - 1)) - 0.5 * (self.rho_s - 1)) - self.omega_c

    @property
    def gamma_b(self):
        """
        UHX 13.5.5
        :return: gamma_b
        """
        if self.config == cfg.a:
            return 0
        elif self.config == cfg.b:
            return (self.G_c - self.C) / self.D_o
        elif self.config == cfg.c:
            return (self.G_c - self.G_1) / self.D_o
        else: # config == d
            return (self.G_c - self.G_s) / self.D_o

    def __step_5_results(self):
        return {
            "gamma": self.gamma,
            "omega_s": self.omega_s,
            "omega_s_star": self.omega_s_star,
            "omega_c": self.omega_c,
            "omega_c_star": self.omega_c_star,
            "gamma_b": self.gamma_b
        }

    @property
    def P_s_prime(self):
        """
        UHX 13.5.6
        :return: Effective shell side pressure
        """
        factor = self.x_s
        factor += 2 * (1 - self.x_t) * self.vu_t
        factor += ((2 * self.vu_s) / self.K_s_t) * (self.D_s / self.D_o) ** 2
        factor -= ((self.rho_s ** 2) - 1) / (self.J * self.K_s_t)
        if self.has_expansion_joint:
            factor -= ((1 - self.J) * (self.D_J ** 2) - self.D_s ** 2) / (2 * self.J * self.K_s_t * self.D_o ** 2)
        return factor * self.P_s

    @property
    def P_t_prime(self):
        """
         UHX 13.5.6
        :return: Effective tube side pressure
        """
        factor = self.x_t
        factor += 2 * (1 - self.x_t) * self.vu_t
        factor += 1 / (self.J * self.K_s_t)
        return factor * self.P_t

    @property
    def P_gamma(self):
        """
         UHX 13.5.6
        :return: Effective pressure due to differential thermal expansion
        """
        return ((self.N_t * self.K_t) / (math.pi * self.a_o ** 2)) * self.gamma

    @property
    def P_w(self):
        """
         UHX 13.5.6
        :return: effective pressure of the mating flange load
        """
        return -1 * (self.U / self.a_o ** 2) * (self.gamma_b / (2 * math.pi)) * self.W_star

    @property
    def P_rim(self):
        """
         UHX 13.5.6
        :return: Effective tubesheet pressure at the rim
        """
        return -1 * (self.U / self.a_o ** 2) * (self.omega_s * self.P_s - self.omega_c * self.P_t)

    @property
    def P_e(self):
        """
         UHX 13.5.6
        :return: Effective pressure on the tubesheet
        """
        a = self.J * self.K_s_t
        b = self.Q_Z_1 + (self.rho_s - 1) * self.Q_Z_2
        c = self.P_s_prime - self.P_t_prime + self.P_gamma + self.P_w + self.P_rim
        return a / (1 + a * b * c)

    def __step_6_results(self):
        return {
            "P_s_prime": self.P_s_prime,
            "P_t_prime": self.P_t_prime,
            "P_gamma": self.P_gamma,
            "P_w": self.P_w,
            "P_rim": self.P_rim,
            "P_e": self.P_e
        }

    @property
    def Q_2(self):
        """
        UHX 13.5.7
        :return: Q_2
        """
        return ((self.omega_s_star * self.P_s - self.omega_c_star * self.P_t) + ((self.gamma_b * self.W_star) / (2 * math.pi))) / (1 + self.phi * self.Z_m)

    @property
    def Q_3(self):
        """
        UHX 13.5.7
        :return:
        """
        if self.P_e == 0:
            RuntimeWarning("Accessed the Q_3 property when P_e was 0")
            return None
        else:
            return self.Q_1 + (2 * self.Q_2 / (self.P_e * self.a_o ** 2))

    def __get_Fm_Step_7(self):
        table = Table_13_1(self.X_a, self.vu_star, self.Q_3)
        self.__F_m = float(table.F_max)

    @property
    def F_m(self):
        """
        UHX 13.5.7
        :return: F_max from table UHX 13.1
        """
        if self.__F_m is None:
            self.__get_Fm_Step_7()

        return self.__F_m

    @property
    def sigma(self):
        """
        UHX 13.5.7
        :return: Max tubesheet bending stress
        """
        if self.P_e == 0:
            return (6 * self.Q_2) / (self.mu_star * ((self.h - self.h_prime_g) ** 2))
        else:
            return (1.5 * self.F_m / self.mu_star) * (((2 * self.a_o) / (self.h - self.h_prime_g)) ** 2) * self.P_e

    @property
    def sigma_acceptable(self):
        """
        UHX 13.5.7
        :return: Acceptibility of tubesheet max bending stress
        """
        if self.is_operating_case:
            return math.fabs(self.sigma) <= self.S_PS
        else:
            return math.fabs(self.sigma) <= 1.5 * self.S

    def __step_7_results(self):
        return {
            "Q_2": self.Q_2,
            "Q_3": self.Q_3,
            "sigma": self.sigma,
            "sigma_acceptable": self.sigma_acceptable
        }

    @property
    def tau_needed(self):
        """
        UHX 13.5.8
        :return: Checks if shear stress check is required
        """
        cond = abs(self.P_e)
        test = (1.6 * self.S * self.mu * self.h) / self.a_o
        return cond > test

    @property
    def tau(self):
        """
        UHX 13.5.8
        :return: Shear stress
        """
        if self.tau_needed:
            return (1 / (4 * self.mu))*((1/ self.h) * (4 * self.A_p / self.C_p))*self.P_e
        else:
            return None

    @property
    def tau_acceptable(self):
        """
        UHX 13.5.8
        :return: determines the acceptibility of shear stress in step 8
        """
        if self.tau_needed:
            return math.fabs(self.tau) <= 0.8 * self.S
        else:
            return None

    def __step_8_results(self):
        return {
            "tau_needed": self.tau_needed,
            "tau": self.tau,
            "tau_acceptable": self.tau_acceptable
        }

    def __get_Ftmin_and_Ftmax(self):
        table = Table_13_2(self.X_a, self.vu_star, self.Q_3, self.P_e)
        self.__F_t_max = float(table.F_t_max)
        self.__F_t_min = float(table.F_t_min)

    @property
    def F_t_max(self):
        """
        UHX 13.5.9
        :return: F_t_max
        """
        if not self.__F_t_max:
            self.__get_Ftmin_and_Ftmax()
        return self.__F_t_max

    @property
    def F_t_min(self):
        """
        UHX 13.5.9
        :return: F_t_min
        """
        if not self.__F_t_min:
            self.__get_Ftmin_and_Ftmax()
        return self.__F_t_min

    @property
    def sigma_t_1(self):
        """
        UHX 13.5.9(a)(1)
        :return: sigma_t_1
        """
        term1 = (1 / (self.x_t - self.x_s))
        term2 = ((self.P_s * self.x_s) - (self.P_t * self.x_t))
        if self.P_e == 0:
            return term1 * (term2 - (2 * self.Q_2 * self.F_t_min / self.a_o ** 2))
        else:
            return term1 * (term2 - (self.P_e * self.F_t_min))

    @property
    def sigma_t_2(self):
        """
        UHX 13.5.9(a)(1)
        :return: sigma_t_2
        """
        term1 = (1 / (self.x_t - self.x_s))
        term2 = ((self.P_s * self.x_s) - (self.P_t * self.x_t))
        if self.P_e == 0:
            return term1 * (term2 - (2 * self.Q_2 * self.F_t_max / self.a_o ** 2))
        else:
            return term1 * (term2 - self.P_e * self.F_t_max)

    @property
    def sigma_t_max(self):
        """
        UHX 13.5.9(a)(2)
        :return: sigma_t_max
        """
        return max(abs(self.sigma_t_1), abs(self.sigma_t_2))

    @property
    def sigma_t_max_acceptable(self):
        """
        UHX 13.5.9(a)(2)
        :return: sigma_t_acceptable
        """
        if self.is_operating_case:
            return self.sigma_t_max < self.S_t * 2
        else:
            return self.sigma_t_max < self.S_t

    @property
    def W_t(self):
        """
        UHX 13.5.9(b)
        :return: tube to tubesheet joint load
        """
        return self.sigma_t_max * math.pi * (self.d_t - self.t_t) * self.t_t

    @property
    def W_t_acceptable(self):
        """
        UHX 13.5.9(a)(2)
        :return: W_t_acceptability
        """
        return self.W_t <= self.L

    @property
    def l_t(self):
        """
        UHX 13.5.9(c)(1)
        :return: tube buckling length
        """

        return self.l * self.k

    @property
    def r_t(self):
        """
        UHX 13.5.9(c)(2)
        :return:
        """
        return 0.25 * ((self.d_t ** 2) + (self.d_t - 2 * self.t_t) ** 2) ** 0.25

    @property
    def F_t(self):
        """
        UHX 13.5.9(c)(2)
        :return: F_t
        """
        return self.l_t / self.r_t

    @property
    def C_t(self):
        """
        UHX 13.5.9(c)(2)
        :return: C_t
        """
        return ((2 * self.E_t * math.pi ** 2) / self.S_y_t) ** 0.5

    @property
    def F_s(self):
        """
        UHX 13.5.9(c)(3)
        :return: F_s
        """
        if self.P_e == 0:
            return 1.25
        else:
            return max(3.25 - 0.25 * (self.Z_d + self.Q_3 * self.Z_w) * self.X_a ** 4, 1.25)

    @property
    def S_tb(self):
        """
        UHX 13.5.9(c)(4)
        :return: S_tb
        """
        if self.C_t <= self.F_t:
            return min((self.E_t * math.pi ** 2) / (self.F_s * self.F_t ** 2), self.S_t)
        else:
            return min((self.S_y_t / self.F_s) * (1 - (self.F_t / (2 * self.C_t))), self.S_t)

    @property
    def sigma_t_min(self):
        """
        UHX 13.5.9(c)(5)
        :return: min sigma_t
        """
        return min(self.sigma_t_1, self.sigma_t_2)

    @property
    def sigma_t_min_acceptable(self):
        """
        UHX 13.5.9(c)(5)
        :return: tube stress acceptable
        """
        return abs(self.sigma_t_min) <= self.S_tb

    def __step_9_results(self):
        common_cases = {"F_t_min": self.F_t_min,
            "F_t_max": self.F_t_max,
            "sigma_t_1": self.sigma_t_1,
            "sigma_t_2": self.sigma_t_2,
            "sigma_t_max": self.sigma_t_max,
            "sigma_t_max_acceptable": self.sigma_t_max_acceptable,
            "W_t": self.W_t,
            "W_t_acceptable": self.W_t_acceptable}

        if self.sigma_t_1 < 0 or self.sigma_t_2 < 0:
            compressive_case = {
                "l_t": self.l_t,
                "r_t": self.r_t,
                "F_t": self.F_t,
                "C_t": self.C_t,
                "F_s": self.F_s,
                "S_tb": self.S_tb,
                "sigma_t_min": self.sigma_t_min,
                "sigma_t_min_acceptable": self.sigma_t_min_acceptable
            }
            return {**common_cases, **compressive_case}
        else:
            return common_cases

    @property
    def sigma_s_m(self):
        """
        UHX 13.5.10(a)
        :return: axial shell membrane stress
        """
        term1 = (self.a_o ** 2) /(self.t_s * (self.D_s + self.t_s))
        term2 = self.P_e + ((self.rho_s ** 2) - 1) * (self.P_s - self.P_t)
        term3 = (self.a_s ** 2) / (self.t_s * (self.D_s + self.t_s))

        return term1 * term2 + term3 * self.P_t

    @property
    def sigma_s_m_acceptable(self):
        """
        UHX 13.5.10(a)(2)
        :return: acceptibility of sigma_s_m
        """

        if self.is_operating_case:
            design_ok = (abs(self.sigma_s_m) <= self.S_PS_s)
        else:
            design_ok = (abs(self.sigma_s_m) <= self.S_t * self.E_s_w)

        if not design_ok or self.sigma_s_m > 0:
            return design_ok
        else:
            return abs(self.sigma_s_m) < self.S_s_b

    def __step_10_results(self):
        return {
            "sigma_s_m": self.sigma_s_m,
            "sigma_s_m_acceptable": self.sigma_s_m_acceptable
        }

    @property
    def L_s_min(self):
        """
        UHX 13.5.11 (a)
        :return: Min length of shell adjacent to tubesheet
        """
        return 1.8 * (self.D_s * self.t_s) ** 0.5

    @property
    def sigma_s_b(self):
        """
        UHX 13.5.11 (a)
        :return: axial bending stress in the shell
        """
        term1 = (6 * self.k_s / self.t_s ** 2)
        term2 = self.beta_s * self.delta_s * self.P_s
        term3 = (6 * (1 - self.vu_star ** 2)) / self.E_star
        term4 = (self.a_o ** 3) / (self.h ** 3)
        term5 = 1 + 0.5 * self.h * self.beta_s
        term6 = self.P_e * (self.Z_v + self.Z_m * self.Q_1)
        term7 = (2 * self.Z_m * self.Q_2) / (self.a_o ** 2)
        return term1 * (term2 + term3 * term4 * term5 * (term6 + term7))

    @property
    def sigma_s(self):
        """
        UHX 13.5.11 (a)
        :return: total shell stress
        """
        return abs(self.sigma_s_m) + abs(self.sigma_s_b)

    @property
    def L_c_min(self):
        """
        UHX 13.5.11 (b)
        :return: Min length of channel adjacent to the tubesheet
        """
        if self.config != cfg.a:
            return 0
        else:
            return 1.8 * (self.D_c * self.t_c) * 0.5

    @property
    def sigma_c_m(self):
        """
        UHX 13.5.11 (b)
        :return: membrane stress in the channel
        """
        return ((self.a_o ** 2) / (self.t_c * (self.D_c + self.t_c))) * self.P_t

    @property
    def sigma_c_b(self):
        """
        UHX 13.5.11 (b)
        :return: bending stress in channel
        """
        term1 = (6 * self.k_c / self.t_c ** 2)
        term2 = self.beta_c * self.delta_c * self.P_t
        term3 = (6 * (1 - self.vu_star ** 2)) / self.E_star
        term4 = (self.a_o ** 3) / (self.h ** 3)
        term5 = 1 + 0.5 * self.h * self.beta_c
        term6 = self.P_e * (self.Z_v + self.Z_m * self.Q_1)
        term7 = (2 * self.Z_m * self.Q_2) / (self.a_o ** 2)
        return term1 * (term2 - term3 * term4 * term5 * (term6 + term7))

    @property
    def sigma_c(self):
        """
         UHX 13.5.11 (b)
        :return: max stress in channel
        """
        return abs(self.sigma_c_m) + abs(self.sigma_c_b)



    @property
    def sigma_s_acceptable(self):
        """
        UHX 13.5.11 (c)
        :return: Determines acceptibility of sigma_s
        """
        if not self.is_operating_case:
            return self.sigma_s <= 1.5 * self.S_s
        else:
            return self.sigma_s <= self.S_PS_s

    @property
    def sigma_c_acceptable(self):
        """
         UHX 13.5.11 (c)
        :return: acceptibility of sigma_c
        """
        if not self.is_operating_case:
            return self.sigma_c <= 1.5 * self.S_c
        else:
            return self.sigma_c <= self.S_PS_c

    def __step_11_results(self):
        if self.config == cfg.d: #  Analysis ends at step 10 for config d
            return {}

        shell = {
            "L_s_min": self.L_s_min,
            "sigma_s_b": self.sigma_s_b,
            "sigma_s": self.sigma_s,
            "sigma_s_acceptable": self.sigma_s_acceptable
        }

        if self.config == cfg.a:
            channel = {
                "L_c_min": self.L_c_min,
                "sigma_c_m": self.sigma_c_m,
                "sigma_c_b": self.sigma_c_b,
                "sigma_c": self.sigma_c,
                "sigma_c_acceptable": self.sigma_c_acceptable
            }
            return {**shell, **channel}
        else:
            return shell

    @property
    def results(self):
        inputs = {}
        outputs = {
            1: self.__step_1_results(),
            2: self.__step_2_results(),
            3: self.__step_3_results(),
            4: self.__step_4_results(),
            5: self.__step_5_results(),
            6: self.__step_6_results(),
            7: self.__step_7_results(),
            8: self.__step_8_results(),
            9: self.__step_9_results(),
            10: self.__step_10_results(),
            11: self.__step_11_results()
        }

        return {
            "inputs": inputs,
            "outputs": outputs
        }


class Bundle:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def add(self, **kwargs):
        self.__dict__.update(kwargs)

if __name__ == "__main__":
    import UHX_11
    p = Bundle(
        r_o=11.5,
        d_t=1,
        t_t=0.049,
        E_t_T=26200000,
        E=26200000,
        S_t_T=23200,
        S=24500,
        p=1.25,
        A_L=0,
        h_g=0,
        CA_t=0,
        h=2)
    p.add(l_x = (0.95 * p.h))
    p.add(pitch_type = UHX_11.PitchType.TRIANGLE)

    tlo_calc_pasrams = UHX_11.UHX11Params(
        radius_to_outermost_tube_hole_center= p.r_o,
        nominal_tube_OD=p.d_t,
        nominal_tube_wall_thickness=p.t_t,
        modulus_of_elasticity_of_tubes_at_tubesheet_design_temperature=p.E_t_T,
        modulus_of_elasticity_for_tubesheet_material_at_tubesheet_design_temperature=p.S,
        allowable_stress_of_tubes_at_tubesheet_design_temperature=p.S_t_T,
        allowable_stress_for_tubesheet_material_at_tubesheet_design_temperature=p.S,
        tube_pitch=p.p,
        area_of_untubed_lanes=p.A_L,
        tube_side_pass_partition_groove_depth=p.h_g,
        tubesheet_corrosion_tube_side=p.CA_t,
        tubesheet_thickness=p.h,
        expanded_depth_of_tube_in_tubesheet=p.l_x,
        pitch_type=p.pitch_type
    )

    tlo_calcs = UHX_11.UHX11(tlo_calc_pasrams)
    p.add(D_o=tlo_calcs.D_o())
    p.add(mu=tlo_calcs.mu())
    p.add(mu_star=tlo_calcs.mu_star())
    p.add(h_prime_g=tlo_calcs.hg_prime())
    p.add(vu_star=tlo_calcs.vu_star())
    p.add(E_star=tlo_calcs.E_star())
    p.add(A=26)
    p.add(C=0)
    p.add(S_PS=0)
    p.add(W_star=0)
    p.add(config=cfg.a)
    p.add(is_operating_case=False)
    p.add(P_s=235)
    p.add(P_t=260)
    p.add(T_t_m=150)
    p.add(E_t=26200000)
    p.add(T_a=70)
    p.add(T_s_m=375)
    p.add(G_c=0)
    p.add(G_s=0)
    p.add(G_1=0)
    p.add(D_s=24)
    p.add(t_s=0.375)
    p.add(E_s=26900000)
    p.add(L=236)
    p.add(vu_s=0.3)
    p.add(alpha_s_m=0.00000705)
    p.add(S_PS_s=0)
    p.add(E_s_w=0)
    p.add(S_s_b=0)
    p.add(S_s=19700)
    p.add(D_c=24)
    p.add(vu_c=0.3)
    p.add(t_c=1)
    p.add(E_c=27900000)
    p.add(S_c=20000)
    p.add(S_PS_c=0)
    p.add(N_t=240)
    p.add(alpha_t_m=00.0000072)
    p.add(vu_t=0.3)
    p.add(A_p=0)
    p.add(C_p=0)
    p.add(L_max=20000)
    p.add(S_t=23200)
    p.add(k=.6)
    p.add(l=26.7647)
    p.add(S_y_t=48800)
    p.add(has_expansion_joint=False)
    p.add(K_J=6000000)


    hx_calcs = UHX_13_5_calcs(p)
    res = hx_calcs.results


    def pretty(d, indent=0):
        for key, value in d.items():
            print('\t' * indent + str(key))
            if isinstance(value, dict):
                pretty(value, indent + 1)
            else:
                print('\t' * (indent + 1) + str(value))
    pretty(res)