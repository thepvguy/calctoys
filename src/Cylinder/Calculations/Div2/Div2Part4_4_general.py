from . import Div2Annex3_D


# Safety Factor, 4.4.2

# Para 4.4.2 equations 4.4.1-3
def FS(F_ic_, S_y):
    if F_ic_ <= 0.55 * S_y:
        return 2
    elif F_ic_ < S_y:
        return 2.407 - 0.714 * (F_ic_ / S_y)
    elif F_ic_ == S_y:
        return 1.667
    else:
        return 0


# Material properties, 4.4.3

# Para 4.4.3.2-1 equation 4.4.4
# elastic ratio
def A_e(F_e, E):
    return F_e / E


def F_ic(F_e, E, sigma_ys, sigma_uts, material_type):
    A_e_ = A_e(F_e, E)
    upper = E
    lower = 0.0
    sigma_guess = 0.0
    tolerance = 0.000001
    diff = 1.0
    while diff > tolerance:
        sigma_guess = 0.5 * (upper + lower)
        tangent_modulus_guess = Div2Annex3_D.E_t(sigma_guess, E, sigma_ys, sigma_uts, material_type)
        inelastic_ratio = sigma_guess / tangent_modulus_guess
        diff = inelastic_ratio - A_e_
        if diff < 0:
            lower = sigma_guess
        else:
            upper = sigma_guess
        diff = abs(diff)

    return min(sigma_guess,
               sigma_ys)  # F_ic is the inelastic buckling stress. Once the predicted inelastic stress is greater than the yield stress, it's no longer inelastic.
    # return sigma_guess
