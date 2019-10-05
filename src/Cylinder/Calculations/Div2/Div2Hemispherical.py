from . import Div2Part4_4_general


# External Pressure, 4.4.7.1


# Para 4.4.7.1-2 equation 4.4.48
def F_he(E_y, t, R_o):
    return 0.075 * E_y * (t / R_o)


# Para 4.4.7.1-2 equation 4.4.50
def F_ha(E_y, R_o, t, sigma_ys, sigma_uts, material_type):
    F_he_ = F_he(E_y, t, R_o)
    F_ic_ = Div2Part4_4_general.F_ic(F_he_, E_y, sigma_ys, sigma_uts, material_type)
    FS_ = Div2Part4_4_general.FS(F_ic_, sigma_ys)

    return F_ic_ / FS_


# Para 4.4.7.1-2 equation 4.4.49
def P_a(E_y, R_o, t, sigma_ys, sigma_uts, material_type):
    F_ha_ = F_ha(E_y, R_o, t, sigma_ys, sigma_uts, material_type)

    return 2 * F_ha_ * (t / R_o)
