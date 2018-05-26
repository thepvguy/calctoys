import math

from ._UHX_common import Configuration as Configuration
from ._UHX_common import AttachmentType as AttachmentType
from .Table_13_1_and_2 import Table_13_1, Table_13_2

class UHX_13_params:
    def __init__(self):
        pass

class UHX_13_calcs:
    def __init__(self, params):

        self.UHX11 = params.UHX11

        self.D_o = self.UHX11.D_o()
        self.isOperatingLoadCase = params.isOperatingLoadCase

        if self.isOperatingLoadCase:
            self.hg_prime = 0
        else:
            self.hg_prime = self.UHX11.hg_prime()
        self.a_s = params.a_s
        self.a_c = params.a_c
        self.N_t = params.N_t
        self.d_t = params.d_t
        self.t_t = params.t_t
        self.t_s = params.t_s
        self.D_s = params.D_s
        self.E_s = params.E_s
        self.E_t = params.E_t
        self.L = params.L

        self.hasExpansionJoint = params.hasExpansionJoint
        self.K_J = params.K_J
        self.vu_s = params.vu_s
        self.vu_c = params.vu_c
        self.vu_t = params.vu_t
        self.config = params.config
        self.channelAttachment = params.channelAttachment
        self.h = params.h
        self.D_c = params.D_c
        self.t_c = params.t_c
        self.E_c = params.E_c
        self.vu_star = self.UHX11.vu_star()
        self.E_star = self.UHX11.E_star()
        self.A = params.A
        self.E = params.E
        self.a_tm = params.a_tm
        self.T_tm = params.T_tm
        self.T_a = params.T_a
        self.a_sm = params.a_sm
        self.T_sm = params.T_sm
        self.G_c = params.G_c
        self.C = params.C
        self.G_1 = params.G_1
        self.G_s = params.G_s
        self.P_s = params.P_s
        self.P_t = params.P_t
        self.W_star = params.W_star
        self.D_J = params.D_J
        self.mu_star = self.UHX11.mu_star()
        self.S = params.S
        self.S_ps = params.S_ps
        self.mu = self.UHX11.mu()
        self.A_p = params.A_p
        self.C_p = params.C_p
        self.S_t = params.S_t
        self.k = params.k
        self.l = params.l
        self.S_y_t = params.S_y_t

    def a_o(self):
        """

        :return: Equivalent radius of OTL
        """
        return self.D_o / 2.0

    def rho_s(self):
        """

        :return: Ratio of shell ID to OTL radius
        """
        return self.a_s / self.a_o()

    def rho_c(self):
        """

        :return: Ratio of channel ID to OTL radius
        """
        return self.a_c / self.a_o

    def x_s(self):
        """

        :return: Something to do with stiffness
        """
        return 1 - self.N_t * (self.d_t / (2.0 * self.a_o())) ** 2

    def x_t(self):
        """

        :return: Something to do with stiffness
        """
        return 1.0 - self.N_t * ((self.d_t - 2 * self.t_t) / (2 * self.a_o())) ** 2


    def K_s(self):
        """

        :return: Shell spring coefficient
        """
        return (math.pi * self.t_s * (self.D_s + self.t_s) * self.E_s) / self.L

    def K_t(self):
        """

        :return: Single tube spring coefficient
        """
        return (math.pi * self.t_t * (self.d_t - self.t_t) * self.E_t) / self.L

    def K_s_t(self):
        """

        :return: Stiffness ratio of bundle to shell
        """
        return self.K_s() / (self.N_t * self.K_t())

    def J(self):
        """

        :return: Expansion joint to shell stiffness ratio
        """
        if self.hasExpansionJoint:
            J = 1 / (1 + (self.K_s() / self.K_J))
        else:
            J = 1
        return J

    def beta_s(self):
        """

        :return: Shell coefficient
        """

        if self.config in [Configuration.a, Configuration.b, Configuration.c]:
            num = (12.0 * (1.0 - self.vu_s ** 2)) ** 0.25
            den = math.sqrt((self.D_s + self.t_s) * self.t_s)
            beta_s = num / den
        elif self.config == Configuration.d:
            beta_s = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return beta_s

    def k_s(self):
        """

        :return: Shell coefficient
        """
        if self.config in [Configuration.a, Configuration.b, Configuration.c]:
            num = self.beta_s() * self.E_s * self.t_s ** 3
            den = 6 * (1 - self.vu_s ** 2)
            k_s = num / den
        elif self.config == Configuration.d:
            k_s = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return k_s

    def lambda_s(self):
        """

        :return: Shell coefficient
        """
        if self.config in [Configuration.a, Configuration.b, Configuration.c]:

            lambda_s = (6 * self.D_s / self. h ** 3) * self.k_s() * (1 + self.h * self.beta_s() + 0.5* ((self.h ** 2) * (self.beta_s() ** 2)))
        elif self.config == Configuration.d:
            lambda_s = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return lambda_s

    def delta_s(self):
        """

        :return: Shell coefficient
        """
        if self.config in [Configuration.a, Configuration.b, Configuration.c]:

            delta_s = ((self.D_s ** 2) / (4 * self.E_s * self.t_s)) * (1 - 0.5 * self.vu_s)
        elif self.config == Configuration.d:
            delta_s = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return delta_s

    def beta_c(self):
        """

        :return: Channel coefficient
        """

        if self.config in [Configuration.b, Configuration.c, Configuration.d]:
            beta_c = 0
        elif self.config == Configuration.a:
            num = (12.0 * (1.0 - self.vu_c ** 2)) ** 0.25
            den = math.sqrt((self.D_c + self.t_c) * self.t_c)
            beta_c = num / den
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return beta_c

    def k_c(self):
        """

        :return: Channel coefficient
        """
        if self.config in [Configuration.b, Configuration.c, Configuration.d]:
            k_c = 0
        elif self.config == Configuration.a:
            num = self.beta_c() * self.E_c * self.t_c ** 3
            den = 6 * (1 - self.vu_c ** 2)
            k_c = num / den
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return k_c

    def lambda_c(self):
        """

        :return: Channel coefficient
        """
        if self.config in [Configuration.b, Configuration.c, Configuration.d]:

            lambda_c = (6 * self.D_c / self. h ** 3) * self.k_c() * (1 + self.h * self.beta_c() + 0.5* ((self.h ** 2) * (self.beta_c() ** 2)))
        elif self.config == Configuration.d:
            lambda_c = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return lambda_c

    def delta_c(self):
        """

        :return: Shell coefficient
        """
        if self.config in [Configuration.a, Configuration.b, Configuration.c]:
            if self.channelAttachment == AttachmentType.CYLINDER:
                delta_c = ((self.D_c ** 2) / (4 * self.E_c * self.t_c)) * (1 - 0.5 * self.vu_c)
            elif self.channelAttachment == AttachmentType.HEMI_HEAD:
                delta_c = ((self.D_c ** 2) / (4 * self.E_c * self.t_c)) * (0.5 * (1 - self.vu_c))
            else:
                raise ValueError("Invalid value for param 'ChannelAttachmentType' <<%r>>" % self.channelAttachment)

        elif self.config == Configuration.d:
            delta_c = 0
        else:
            raise ValueError("Invalid value for param 'Configuration' <<%r>>" % self.config)

        return delta_c

    def X_a(self):
        return (
                   24 *
                   (1 - self.vu_star ** 2) *
                   self.N_t *
                    (
                        (self.E_t * self.t_t * (self.d_t - self.t_t) * self.a_o() ** 2) /
                        (self.E_star * self.L * self. h ** 3)
                    )
               ) ** 0.25

    def K(self):
        return self.A / self.D_o

    def F(self):
        return ((1 - self.vu_star) / self.E_star) * (self.lambda_s() + self.lambda_c() + self.E * math.log(self.K()))

    def __table13_1_(self):
        # CAUTION! q_3 is NOT set, trying to access F_m before setting it will cause a crash!
        return Table_13_1(self.X_a(), self.vu_star)

    def phi(self):
        return (1 + self.vu_star) * self.F()

    def Q_1(self):
        table = self.__table13_1_()
        return (self.rho_s() - 1 - self.phi() * table.Z_v) / (1 + self.phi() * table.Z_m)

    def Q_Z1(self):
        table = self.__table13_1_()
        return 0.5 * (table.Z_d + self.Q_1() * table.Z_m * self.X_a() ** 4)

    def Q_Z2(self):
        table = self.__table13_1_()
        return 0.5 * (table.Z_v + self.Q_1() * table.Z_m) * self.X_a() ** 4

    def U(self):
        table = self.__table13_1_()
        return ((table.Z_w + (self.rho_s() - 1) * table.Z_m) * self.X_a() ** 4) / (1 + self.phi() * table.Z_m)

    def gamma(self):
        if self.isOperatingLoadCase:
            gamma = (self.a_tm * (self.T_tm - self.T_a) - self.a_sm * (self.T_sm - self.T_a)) * self.L
        elif not self.isOperatingLoadCase:
            gamma = 0
        else:
            raise ValueError("Parameter isOperatingLoadCase is not set!")
        return gamma

    def omega_s(self):
        return self.rho_s() * self.k_s() * self.beta_s() * self.delta_s() * (1 - self.h * self.beta_s())

    def omega_s_star(self):
        return (self.a_o() ** 2) * 0.25 * ((self.rho_s() ** 2) - 1) * (self.rho_s() -1 ) - self.omega_s()

    def omega_c(self):
        return self.rho_c() * self.k_c() * self.beta_c() * self.delta_c() * (1 + self.h * self.beta_c())

    def omega_c_star(self):
        return (self.a_o() ** 2) *(0.25 * ((self.rho_c() ** 2) - 1) * (self.rho_c() -1 ) - 0.5 * (self.rho_s() - 1)) - self.omega_c()

    def gamma_b(self):
        if self.config == Configuration.a:
            gamma_b = 0
        elif self.config == Configuration.b:
            gamma_b = (self.G_c - self.C) / self.D_o
        elif self.config == Configuration.c:
            gamma_b = (self.G_c - self.G_1) / self.D_o
        elif self.config == Configuration.d:
            return (self.G_c - self.G_s) /self.D_o
        else:
            raise ValueError("Invalid configuration <<%r>> for UHX 13 HX" % self.config)

        return gamma_b

    def P_s_prime(self):
        """

        :return: Effective Shell Side Pressure
        """
        term1 = self.x_s()
        term2 = 2 * (1- term1) * self.vu_t
        term3 = (2 / self.K_s_t()) * ((self.D_s / self.D_o) ** 2) * self.vu_s
        term4 = -((self.rho_s() ** 2) - 1) / (self.J() * self.K_s_t())
        term5 = -(1 - self.J()) / (2 * self.K_s_t())
        term6 = ((self.D_J ** 2)  - self.D_s ** 2) / (self.D_o ** 2)

        return self.P_s * (term1 + term2 + term3 + term4 + term5 * term6)

    def P_t_prime(self):
        """

        :return: Effective tube side pressure
        """
        return self.P_t * (self.x_t() + 2 * (1 - self.x_t()) * self.vu_t + (1 / (self.J() * self.K_s_t())))

    def P_gamma(self):
        """

        :return: Effective tube pressure
        """
        return self.gamma() * ((self.N_t * self.K_t()) / (math.pi * (self.a_o() ** 2)))

    def P_w(self):
        """

        :return: Effective pressure load from mating component
        """
        return - self.W_star * ((self.U() * self.gamma_b()) / ((self.a_o() ** 2) * 2 * math.pi))

    def P_rim(self):
        return  - (self.U() / self.a_o() ** 2) * (self.omega_s_star() * self.P_s - self.omega_c_star() * self.P_t)

    def P_e(self):
        term1 = self.J() * self.K_s_t()
        term2 = self.J() * self.K_s_t() * (self.Q_Z1() + (self.rho_s() - 1) * self.Q_Z2())
        term3 = (self.P_s_prime() - self.P_t_prime() + self.P_gamma() + self.P_w() + self.P_rim())
        return (term1 / (1 + term2)) * term3

    def Q_2(self):

        return ((self.omega_s_star() * self.P_s - self.omega_c_star() * self.P_t) + (((self.gamma_b())/(2 * math.pi)) * self.W_star)) / (1 + self.phi() * self.__table13_1_().Z_m)

    def Q_3(self):
        return self.Q_1() + ((2 * self.Q_2()) / (self.P_e() * self.a_o() ** 2))

    def sigma_bending(self):
        if self.P_e() != 0:
            table = self.__table13_1_()
            table.q_3 = self.Q_3()

            term1 = 1.5 * table.F_max / self.mu_star
            term2 = ((2 * self.a_o()) / (self.h - self.hg_prime)) ** 2

            stress =  term1 * term2 * self.P_e()
        else:
            stress = (6 * self.Q_2()) / (self.mu_star * (self.h - self.hg_prime) ** 2)
        return stress

    def acceptable_bending(self):
        if self.isOperatingLoadCase:
            return self.S_ps
        else:
            return 1.5 * self.S

    def shear_required(self):
        return abs(self.P_e()) <= ((1.6 * self.S * self.mu * self.h) / (self.a_o()))

    def tau(self):
        return (0.25 / self.mu) * ((1 / self.h) * (4 * (self.A_p / self.C_p))) * self.P_e()

    def __table_13_2(self, has_external_pressure):
        return Table_13_2(
            xa=self.X_a(), v_star=self.vu_star, q_3=self.Q_3(), hasexternalpressure=has_external_pressure
        )

    def sigma_t1(self):
        term1 = (1 / (self.x_t() - self.x_s()))
        term2 = self.P_s * self.x_s() - self.P_t * self.x_t()

        if self.P_e() == 0:
            term3 = (2 * self.Q_2() * self.__table_13_2(False).F_t_min) / self.a_o() ** 2
        else:
            term3 = self.P_e() * self.__table_13_2(True).F_t_min

        return term1 * (term2 - term3)

    def sigma_t2(self):
        term1 = (1 / (self.x_t() - self.x_s()))
        term2 = self.P_s * self.x_s() - self.P_t * self.x_t()

        if self.P_e() == 0:
            term3 = (2 * self.Q_2() * self.__table_13_2(False).F_t_max) / self.a_o() ** 2
        else:
            term3 = self.P_e() * self.__table_13_2(True).F_t_max

        return term1 * (term2 - term3)

    def sigma_t_max(self):
        return max(abs(self.sigma_t1()), abs(self.sigma_t2()))

    def acceptable_sigma_t(self):
        if self.isOperatingLoadCase:
            return 2 * self.S_t
        else:
            return self.S_t

    def W_t(self):
        return self.sigma_t_max() * math.pi * (self.d_t - self.t_t) * self.t_t

    def l_t(self):
        return self.k * self.l

    def r_t(self):
        return 0.25 * math.sqrt((self.d_t ** 2) + (self.d_t - 2 * self.t_t) ** 2)

    def F_t(self):
        return self.l_t() / self.r_t()

    def C_t(self):
        return math.sqrt((2 * (math.pi ** 2) * self.E_t) / (self.S_y_t))

    def F_s(self):
        if self.P_e() == 0:
            return 1.25
        else:
            table = self.__table13_1_()
            fs = 3.25 - 0.25 * (table.Z_d + self.Q_3() * table.Z_w) * self.X_a() ** 4
            return max(fs,1.25)

    def S_tb(self):
        if self.C_t() <= self.F_t():
            stress_limit = ((1 / self.F_s()) * (math.pi ** 2) * (self.E_t / self.F_t() ** 2))
        else:
            stress_limit = (self.S_y_t / self.F_s()) * (1 - (self.F_t() / (2 * self.C_t())))
        return min(stress_limit, self.S_t)

    def sigma_t_min(self):
        return min(self.sigma_t1(), self.sigma_t2())

    def acceptable_sigma_tb(self):
        return self.S_tb()

    # TODO: Step 10 and 11; UHX 13.6, 13.7, 13.8, 13.9
