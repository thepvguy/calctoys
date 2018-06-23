import math


class AppendixO:
    def __init__(
            self,
            root_area,
            gasket_od,
            gasket_id,
            nut_factor,
            number_of_bolts,
            max_design_pressure,
            flange_yield_stress_at_assembly,
            flange_yield_at_operation,
            max_allowable_bolt_stress,
            min_allowable_bolt_stress,
            bolt_yield_stress_at_ambient,
            max_bolt_stress_new,
            target_gasket_stress,
            max_gasket_stress,
            min_gasket_seating_stress,
            min_gasket_operating_stress,
            bolt_diameter,
            total_flange_rotation,
            max_flange_rotation,
            min_fraction_of_gasket_load_after_relaxations,
            is_metric
    ):
        self.A_b = root_area
        self.G_od = gasket_od
        self.G_id = gasket_id
        self.K = nut_factor
        self.n_b = number_of_bolts
        self.P_max = max_design_pressure
        self.S_ya = flange_yield_stress_at_assembly
        self.S_yo = flange_yield_at_operation
        self.Sb_max = max_allowable_bolt_stress
        self.Sb_min = min_allowable_bolt_stress
        self.Sy_b_a = bolt_yield_stress_at_ambient
        self.Sf_max = max_bolt_stress_new
        self.S_g_T = target_gasket_stress
        self.S_g_max = max_gasket_stress
        self.S_g_min_S = min_gasket_seating_stress
        self.S_g_min_O = min_gasket_operating_stress
        self.phi_b = bolt_diameter
        self.theta_f_max = total_flange_rotation
        self.theta_g_max = max_flange_rotation
        self.rho_g = min_fraction_of_gasket_load_after_relaxations
        self.is_metric = is_metric
        self.messages = []
        self.check_inputs()

    def check_inputs(self):
        """

        :return:
        """
        self.messages = []
        # Check O-4.1(a)
        if not 0.3 <= self.theta_g_max <= 1:
            self.messages.append("Flange rotation limits are outside of typical ranges")

        # O-4.1(b)
        if not 0.4 <= self.Sy_b_a <= 0.7:
            self.messages.append("Max permissible bolt stress it outside typical ranges")

        # O-4.1(c)
        if not 0.2 <= self.Sb_min <= 0.4:
            self.messages.append("Minimum permissible bolt stress is outside typical ranges")

    @property
    def A_g(self):
        """
        Defined in O-1.3 Definitions
        :return: gasket area
        """
        return 0.25 * math.pi * ((self.G_od ** 2) - (self.G_id ** 2))

    @property
    def Sb_sel_typ(self):
        """
        Equation O-1
        :return: Assembly bolt stress for a typical joint
        """
        return self.S_g_T * (self.A_g / (self.n_b * self.A_b))

    @property
    def T_b(self):
        """
        Defined in O-2M and O-2
        :return: Assembly bolt torque
        """
        divisor = 1000 if self.is_metric else 12
        return self.Sb_sel_typ * self.K * self.A_b * self.phi_b / divisor

    def T_b_modify_table(self, Sb_sel_prime, T_b_prime, K_prime):
        """
        An example of the type of table produced using equations O-2M and O-2
        method is given in Table 1, which was constructed using
        a  bolt  stress  of  50  ksi  and  a  nut  factor, K, of 0.20.

        If another  bolt  stress  or  nut  factor  is  required,  then  the
        table may be converted to the new values using eq. (O-3),
        where Sb_sel_prime, T_b_prime, and K_prime are the original values.
        :return: Modified Table 1 value
        """

        if not (Sb_sel_prime and T_b_prime and K_prime):
            ValueError("Required parameter is zero or none")

        return (self.K * self.Sb_sel_typ) / (K_prime * Sb_sel_prime * T_b_prime)

    @property
    def Sb_sel_upper(self):
        """
        Equations O-4
        :return: Either typical or upper limit on bolt stress
        """
        return min(self.Sb_sel_typ, self.Sb_max)

    @property
    def Sb_sel_lower(self):
        """
        Equations O-5
        :return: Either typical or lower limit on bolt stress
        """
        return max(self.Sb_sel_typ, self.Sb_min)

    @property
    def Sb_sel_flange_limit(self):
        """
        Equations O-6
        :return: Either typical or lower limit on flange stress
        """
        return min(self.Sb_sel_typ, self.Sf_max)

    @property
    def Sb_sel(self):
        return max(self.Sb_sel_upper, self.Sb_sel_lower, self.Sb_sel_flange_limit)

    @property
    def seating_stress_limit(self):
        """
        Right hand side of O-7
        :return:
        """
        return self.S_g_min_S * (self.A_g / (self.A_b * self.n_b))

    @property
    def min_assembly_seating_stress_check(self):
        """
        Equation O-7
        :return:
        """
        return self.Sb_sel >= self.seating_stress_limit

    @property
    def min_operating_seating_stress_limit(self):
        """
        right hand side of O-8
        :return:
        """
        return (self.S_g_min_O * self.A_g + 0.25 * self.P_max * self.G_id ** 2) / (self.rho_g * self.A_b * self.n_b)

    @property
    def min_operating_seating_stress_check(self):
        """
        Equation O-8
        :return:
        """
        return self.Sb_sel >= self.min_operating_seating_stress_limit

    @property
    def gasket_max_stress_limit(self):
        """
        right hand side of O-9
        :return:
        """
        return self.S_g_max * (self.A_g / (self.A_b * self.n_b))

    @property
    def max_gasket_stress_check(self):
        """
        equation O-9
        :return:
        """
        return self.Sb_sel <= self.gasket_max_stress_limit

    @property
    def flange_rotation_limit(self):
        """
        right hand side of O-10
        :return:
        """
        return self.Sf_max * (self.theta_g_max / self.theta_f_max)

    @property
    def max_flange_rotations_check(self):
        """
        equation O-10
        :return:
        """
        return self.Sb_sel <= self.flange_rotation_limit

    def get_result_dict(self):
        return {
            "input": {
                "A_b": self.A_b,
                "G_od": self.G_od,
                "G_id": self.G_id,
                "K": self.K,
                "n_b": self.n_b,
                "P_max": self.P_max,
                "S_ya": self.S_ya,
                "S_yo": self.S_yo,
                "Sb_max": self.Sb_max,
                "Sb_min": self.Sb_min,
                "Sy_b_a": self.Sy_b_a,
                "Sf_max": self.Sf_max,
                "S_g_T": self.S_g_T,
                "S_g_max": self.S_g_max,
                "S_g_min_S": self.S_g_min_S,
                "S_g_min_O": self.S_g_min_O,
                "phi_b": self.phi_b,
                "theta_f_max": self.theta_f_max,
                "theta_g_max": self.theta_g_max,
                "rho_g": self.rho_g,
                "is_metric": self.is_metric,
            },
            "output": {
                "A_g": self.A_g,
                "Sb_sel_typ": self.Sb_sel_typ,
                "T_b": self.T_b,
                "Sb_sel_upper": self.Sb_sel_upper,
                "Sb_sel_lower": self.Sb_sel_lower,
                "Sb_sel_flange_limit": self.Sb_sel_flange_limit,
                "Sb_sel": self.Sb_sel,
                "seating_stress_limit": self.seating_stress_limit,
                "min_operating_seating_stress_limit": self.min_operating_seating_stress_limit,
                "gasket_max_stress_limit": self.gasket_max_stress_limit,
                "flange_rotation_limit": self.flange_rotation_limit
            },
            "acceptance": {
                "min_assembly_seating_stress_check": self.min_assembly_seating_stress_check,
                "min_operating_seating_stress_check": self.min_operating_seating_stress_check,
                "max_gasket_stress_check": self.max_gasket_stress_check,
                "max_flange_rotations_check": self.max_flange_rotations_check
            },
            "messages": self.messages
        }
