import math


class RCB_4:
    def __init__(self, shell_ID, tube_pitch, tube_OD, layout_pattern, fluid_density, flow_rate, max_rhovsquared, hasImpingement):
        self.D_s = shell_ID
        self.P_t = tube_pitch
        self.D_t = tube_OD
        self.pat = layout_pattern
        self.rho = fluid_density
        self.Q = flow_rate
        self.rho_v_squared = max_rhovsquared
        self.hasImpingement = hasImpingement
        self.V = math.sqrt(self.rho_v_squared / self.rho)
        self.A_s = (self.Q / (self.rho * self.V)) * 0.04 # 0.04 is the unit conversion from (ft^2 * s) / hr to in^2; (144 / 3600)

    @property
    def F_1(self):
        if self.hasImpingement == True:
            return 0
        elif self.hasImpingement == False:
            return 1
        else:
            raise (ValueError("Impingement presence not set"))

    @property
    def F_2(self):
        if self.pat in [layout_patterns.square, layout_patterns.rotated_triangle]:
            return 1
        elif self.pat == layout_patterns.triangle:
            return 0.866
        elif self.pat == layout_patterns.rotated_square:
            return 0.707
        else:
            raise ValueError("Invalid tube pattern type")

    def min_free_height_62(self, nozzle_ID):
        print("TEMA RGP RCB-4.62")
        print("(As - F1(π / 4 * Dn^2) * (Pt - Dt) / (F2 * Pt)) / (π * Dn)")
        print("(%f - %f * (pi / 4 * %f^2 ) * (%f - %f) / (%f * %f )) / (pi * %f)" % (self.A_s, self.F_1, nozzle_ID, self.P_t, self.D_t, self.F_2, self.P_t, nozzle_ID))
        return (self.A_s - self.F_1 * (math.pi / 4 * (nozzle_ID ** 2)) * (self.P_t - self.D_t) / (self.F_2 * self.P_t)) / (math.pi * nozzle_ID)


class layout_patterns(enumerate):
    triangle = 30
    square = 90
    rotated_triangle = 60
    rotated_square = 45

if __name__ == "__main__":
    calc = RCB_4(36, .7813, .625, layout_patterns.rotated_square, 64.4, 4000, 4000, False)
    print("Inlet nozzle:")
    print(calc.min_free_height_62(12))
    print("Outlet nozzle:")
    print(calc.min_free_height_62(6.219))