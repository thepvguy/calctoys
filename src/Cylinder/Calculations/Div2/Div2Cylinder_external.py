import math
import enum

from . import CalcExceptions, Div2Part4_4_general


# Values given in nomenclature section
class CompressionMemberRotationalCoefficient(enum.Enum):
    SIDEWAYS_COMPRESSION = enum.auto()
    SHEAR_NOT_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED_SINGLE_CURVATURE = enum.auto()
    SHEAR_NOT_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED_DOUBLE_CURVATURE = enum.auto()
    SHEAR_NOT_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED = enum.auto()
    SHEAR_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED = enum.auto()
    COMPRESSION_ROTATIONALLY_UNBRACED = enum.auto()
    UNBRACED_SKIRT_SUPPORTED = enum.auto()


def C_m(C_m_type, M_1=0, M_2=0):
    if C_m_type == CompressionMemberRotationalCoefficient.SIDEWAYS_COMPRESSION:
        return 0.85
    elif C_m_type == CompressionMemberRotationalCoefficient.SHEAR_NOT_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED_SINGLE_CURVATURE:
        if M_1 > M_2:
            return 0.6 - 0.4 * (-abs(M_1 / M_2))
        else:
            return 0.6 - 0.4 * (-abs(M_2 / M_1))
    elif C_m_type == CompressionMemberRotationalCoefficient.SHEAR_NOT_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED_DOUBLE_CURVATURE:
        if M_1 > M_2:
            return 0.6 - 0.4 * (abs(M_1 / M_2))
        else:
            return 0.6 - 0.4 * (abs(M_2 / M_1))
    elif C_m_type == CompressionMemberRotationalCoefficient.SHEAR_BETWEEN_SUPPORT_ROTATIONALLY_RESTRAINED:
        return 0.85
    elif C_m_type == CompressionMemberRotationalCoefficient.COMPRESSION_ROTATIONALLY_UNBRACED:
        return 1.0
    elif C_m_type == CompressionMemberRotationalCoefficient.UNBRACED_SKIRT_SUPPORTED:
        return 1.0


class AxialCompressionEndCondition(enum.Enum):
    FIXED_FIXED = enum.auto()
    FIXED_FREE = enum.auto()
    PINNED_FIXED = enum.auto()
    PINNED_PINNED = enum.auto()


def K_u(axial_compression_end_condition):
    if axial_compression_end_condition == AxialCompressionEndCondition.FIXED_FIXED:
        return 0.65
    elif axial_compression_end_condition == AxialCompressionEndCondition.FIXED_FREE:
        return 2.10
    elif axial_compression_end_condition == AxialCompressionEndCondition.PINNED_FIXED:
        return 0.80
    elif axial_compression_end_condition == AxialCompressionEndCondition.FIXED_FIXED:
        return 1.00
    else:
        ValueError("Axial compression end condition is not a valid choice")


# External Pressure, 4.4.5.1

# Para 4.4.5.1-2 equation 4.4.18
# Also 4.4.12.2(k) equation 103
def M_x(L, R_o, t):
    return L / math.sqrt(R_o * t)


# Para 4.4.5.1-2 equation 4.4.19 - 4.4.22
def C_h(D_o, t, L):
    M_x_ = M_x(L, D_o / 2, t)

    if M_x_ >= (2 * D_o / t) ** 0.94:
        return 0.55 * (t / D_o)
    elif M_x_ > 13:
        return 1.12 * M_x_ ** -1.058
    elif M_x_ > 1.5:
        return 0.92 / (M_x_ - 0.579)
    else:  # M_x_ <= 1.5
        return 1


# Para 4.4.5.1-1 equation 4.4.17
def F_he(E_y, D_o, t, L):
    C_h_ = C_h(D_o, t, L)

    return (1.6 * C_h_ * E_y * t) / D_o


def FS_ha(D_o, t, L, E_y, S_y, S_u, material_type):
    F_he_ = F_he(E_y, D_o, t, L)
    F_ic_ = Div2Part4_4_general.F_ic(F_he_, E_y, S_y, S_u, material_type)
    return Div2Part4_4_general.FS(F_ic_, S_y)


# Para 4.4.5.1-5 equation 4.4.24
def F_ha(D_o, t, L, E_y, S_y, S_u, material_type):
    F_he_ = F_he(E_y, D_o, t, L)
    F_ic_ = Div2Part4_4_general.F_ic(F_he_, E_y, S_y, S_u, material_type)
    FS_ = Div2Part4_4_general.FS(F_ic_, S_y)

    return F_ic_ / FS_


# Para 4.4.5.1-5 equation 4.4.23
def P_a(E_y, D_o, t, L, S_y, S_u, material_type):
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)

    return 2 * F_ha_ * (t / D_o)


# Section 4.4.12.2, combined loadings and allowable compressive stresses
# 4.4.12.2 (b), Axial Compressive stress acting alone

# 4.4.12.2(b)(1)
def is_local_buckling(lambda_c_):
    return lambda_c_ <= 0.15


# 4.4.12.2(b)(2)
def column_buckling_upper_limit_exceeded(axial_compression_end_condition, L_u, D_o, D_i):
    r_g_ = r_g(D_o, D_i)
    K_u_ = K_u(axial_compression_end_condition)

    return ((K_u_ * L_u) / r_g_) < 200


# 4.4.55-57
def c_bar(L, D_o, t):
    M_x_ = M_x(L, D_o / 2, t)
    if M_x_ >= 15:
        return 1
    elif M_x_ > 1.5:
        return 3.13 / (M_x_ ** 0.42)
    else:
        return 2.64


# 4.4.53-54
def C_x(L, D_o, t):
    c_ = c_bar(L, D_o, t)

    if (D_o / t) < 1247:
        return min(((409 * c_) / (389 + (D_o / t))), 0.9)
    elif (D_o / t) < 2000:
        return 0.25 * c_
    else:
        CalcExceptions.Div2CylinderException(
            f"Value of D_o / t -> {D_o} / {t} -> {D_o / t}, is greater than 2000 and out of range for C_x")


# 4.4.52
def F_xe(E_y, L, D_o, t):
    C_x_ = C_x(L, D_o, t)
    return (C_x_ * E_y * t) / D_o


# 4.4.58
def F_xa(D_o, t, L, E_y, S_y, S_u, material_type):
    F_xe_ = F_xe(E_y, L, D_o, t)
    F_ic_ = Div2Part4_4_general.F_ic(F_xe_, E_y, S_y, S_u, material_type)
    FS_ = Div2Part4_4_general.FS(F_ic_, S_y)

    return F_ic_ / FS_


def FS_xa(D_o, t, L, E_y, S_y, S_u, material_type):
    F_xe_ = F_xe(E_y, L, D_o, t)
    F_ic_ = Div2Part4_4_general.F_ic(F_xe_, E_y, S_y, S_u, material_type)
    return Div2Part4_4_general.FS(F_ic_, S_y)


# 4.4.59 - 60
def F_ca(axial_compression_end_condition, L, D_o, t, E_y, S_y, S_u, material_type, F_xa_=None):
    # part (e) uses this equation with a modified parameter for F_xa
    if F_xa_ is None:
        _F_xa_ = F_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    else:
        _F_xa_ = F_xa_

    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)

    if l_c >= 1.2:
        return 0.88 * _F_xa_ / (l_c ** 2)
    elif l_c > 0.15:
        return _F_xa_ * (1 - 0.74 * (l_c - 0.15)) ** 0.3
    else:
        CalcExceptions.Div2CylinderException("lambda_c is outside range applicable to F_ca")


# 4.4.12.2 (c), Compressive bending stress

# The allowable axial compressive membrane stress
# of a cylindrical shell subject to a bending moment
# acting across the full circular cross section F_ba,
# shall be determined using the procedure in 4.4.12.2(b)

def F_ba(D_o, t, L, E_y, S_y, S_u, material_type):
    return F_xa(D_o, t, L, E_y, S_y, S_u, material_type)


# 4.4.12.2 (d), shear stress
# 4.4.66-67
def alpha_v(D_o, t):
    if (D_o / t) <= 500:
        return 0.8
    else:
        return 1.389 - 0.218 * math.log10(D_o / t)


# 4.4.62-65
def C_v(L, D_o, t):
    M_x_ = M_x(L, D_o / 2, t)

    if M_x_ >= (4.347 * (D_o / t)):
        return 0.716 * math.sqrt(t / D_o)
    elif M_x_ >= 26:
        return 1.492 / math.sqrt(M_x_)
    elif M_x_ > 1.5:
        return (9.64 / M_x_ ** 2) * math.sqrt(1 + 0.0239 * M_x_ ** 3)
    else:
        return 4.454


# 4.4.61
def F_ve(L, D_o, t, E_y):
    a_v = alpha_v(D_o, t)
    C_v_ = C_v(L, D_o, t)

    return a_v * C_v_ * E_y * (t / D_o)


# 4.4.68
def F_va(L, D_o, t, E_y, S_y, S_u, material_type):
    F_ve_ = F_ve(L, D_o, t, E_y)
    F_ic_ = Div2Part4_4_general.F_ic(F_ve_, E_y, S_y, S_u, material_type)
    FS_ = Div2Part4_4_general.FS(F_ic_, S_y)

    return F_ic_ / FS_


# 4.4.12.2 (e), Axial compressive stress and hoop compression
# 4.4.72
def f_x(F, P, D_o, D_i):
    return f_a(F, D_o, D_i) + f_q(P, D_o, D_i)


# 4.4.71
def C_2(F, P, D_o, D_i, l_c=None):
    if l_c is None or l_c <= 0.15:
        num = f_x(F, P, D_o, D_i)
    else:  # l_c > 0.15:
        num = f_a(F, D_o, D_i)
    return num / f_h(P, D_o, (D_o - D_i) / 2)


# 4.4.70
def C_1(D_o, t, L, E_y, S_y, S_u, material_type):
    F_xa_ = F_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    FS_xa_ = FS_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)
    FS_ha_ = FS_ha(D_o, t, L, E_y, S_y, S_u, material_type)

    return ((F_xa_ * FS_xa_ + F_ha_ * FS_ha_) / S_y) - 1


def F_xha_for_e_2(axial_compression_end_condition, P, F, D_o, t, L, E_y, S_y, S_u, material_type):
    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)
    C_1_ = C_1(D_o, t, L, E_y, S_y, S_u, material_type)
    C_2_ = C_2(F, P, D_o, D_o - 2 * t, l_c)
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)
    F_xa_ = F_xa(D_o, t, L, E_y, S_y, S_u, material_type)

    return ((1 / (F_xa_ ** 2)) - (C_1_ / (C_2_ * F_xa_ * F_ha_)) + (1 / (C_2_ ** 2 * F_ha_ ** 2))) ** -0.5


# 4.4.74
def F_ah2(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type):
    F_xa_ = F_xha_for_e_2(axial_compression_end_condition, P, F, D_o, t, L, E_y, S_y, S_u, material_type)
    F_ca_ = F_ca(axial_compression_end_condition, L, D_o, t, E_y, S_y, S_u, material_type, F_xa_)
    f_q_ = f_q(P, D_o, (D_o / 2) - t)

    return F_ca_ * (1 - (f_q_ / S_y))


# 4.4.69 & 73
def F_xha(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type):
    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)

    if l_c >= 1.2:
        CalcExceptions.Div2CylinderException("lambda_c is greater than 1.2, this does not apply")

    F_ah1 = F_xha_for_e_2(axial_compression_end_condition, P, F, D_o, t, L, E_y, S_y, S_u, material_type)

    if l_c <= 0.15:
        return F_ah1

    F_ah2_ = F_ah2(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type)

    return min(F_ah1, F_ah2_)


# 4.4.75
def F_hxa(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type):
    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)
    if l_c <= 0.15:
        F_xha_ = F_xha(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type)
        C_2_ = C_2(F, P, D_o, D_o - 2 * t)

        return F_xha_ / C_2_
    else:
        CalcExceptions.Div2CylinderException("lambda_c is out of range for F_hxa")


# 4.4.12.2 (f),Compressive Bending Stress and hoop compression
# 4.4.79
def n(D_o, t, L, E_y, S_y, S_u, material_type):
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)
    FS_ha_ = FS_ha(D_o, t, L, E_y, S_y, S_u, material_type)

    return 5 - (4 * F_ha_ * FS_ha_) / S_y


# 4.4.77
def C_4(P, M, D_o, t, L, E_y, S_y, S_u, material_type):
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)
    F_ba_ = F_ba(D_o, t, L, E_y, S_y, S_u, material_type)
    f_h_ = f_h(P, D_o, t)
    S_ = S(D_o, D_o - 2 * t)
    f_b_ = f_b(M, S_)

    return f_b_ * F_ha_ / (f_h_ * F_ba_)


# 4.4.78
def C_3(P, M, D_o, t, L, E_y, S_y, S_u, material_type):
    C_4_ = C_4(P, M, D_o, t, L, E_y, S_y, S_u, material_type)
    n_ = n(D_o, t, L, E_y, S_y, S_u, material_type)
    c3_const = C_4_ ** 2 + 0.6 * C_4_

    def newton_raphson_c3(c3):
        # The c3 function
        num = (c3_const * c3 ** 2 + (c3 ** (2 * n_)) - 1)

        # Derivative of the c3 function
        den = 2 * c3_const * c3 + 2 * n_ * c3 ** (2 * n_ - 1)

        return c3 - (num / den)

    # c4 is just a guess
    # We're looking for this to equal zero
    guess_1 = 0
    guess_2 = newton_raphson_c3(C_4_)
    while abs((guess_1 - guess_2)) > 0.0000001:
        guess_1 = guess_2
        guess_2 = newton_raphson_c3(guess_1)

    return guess_2


# 4.4.76
def F_bha(P, M, D_o, t, L, E_y, S_y, S_u, material_type):
    C_3_ = C_3(P, M, D_o, t, L, E_y, S_y, S_u, material_type)
    C_4_ = C_4(P, M, D_o, t, L, E_y, S_y, S_u, material_type)
    F_ba_ = F_ba(D_o, t, L, E_y, S_y, S_u, material_type)

    return C_3_ * C_4_ * F_ba_


# 4.4.80
def F_hba(P, M, D_o, t, L, E_y, S_y, S_u, material_type):
    F_bha_ = F_bha(P, M, D_o, t, L, E_y, S_y, S_u, material_type)
    f_h_ = f_h(P, D_o, t)
    S_ = S(D_o, D_o - 2 * t)
    f_b_ = f_b(M, S_)

    return F_bha_ * (f_h_ / f_b_)


# 4.4.12.2 (g), Shear Stress and Hoop Compression
# 4.4.82
def C_5(V, phi, P, D_o, t):
    f_v_ = f_v(V, phi, D_o, D_o - 2 * t)
    f_h_ = f_h(P, D_o, t)

    return f_v_ / f_h_


# 4.4.81
def F_vha(V, phi, P, L, D_o, t, E_y, S_y, S_u, material_type):
    F_va_ = F_va(L, D_o, t, E_y, S_y, S_u, material_type)
    C_5_ = C_5(V, phi, P, D_o, t)
    F_ha_ = F_ha(D_o, t, L, E_y, S_y, S_u, material_type)

    return ((((F_va_ ** 2) / (2 * C_5_ * F_ha_)) ** 2 + F_va_ ** 2) ** 0.5) - (F_va_ / (2 * C_5_ * F_ha_))


# 4.4.83
def F_hva(V, phi, P, L, D_o, t, E_y, S_y, S_u, material_type):
    C_5_ = C_5(V, phi, P, D_o, t)
    F_vha_ = F_vha(V, phi, P, L, D_o, t, E_y, S_y, S_u, material_type)

    return F_vha_ / C_5_


# 4.4.12.2 (h), Shear Stress and Hoop Compression
# 4.4.84
def K_S(V, phi, L, D_o, t, E_y, S_y, S_u, material_type):
    f_v_ = f_v(V, phi, D_o, D_o - 2 * t)
    F_va_ = F_va(L, D_o, t, E_y, S_y, S_u, material_type)

    return 1.0 - (f_v_ / F_va_) ** 2


# 4.4.89
def F_e(axial_compression_end_condition, D_o, t, L_u, E_y):
    r_g_ = r_g(D_o, D_o - 2 * t)
    K_u_ = K_u(axial_compression_end_condition)

    return ((math.pi ** 2) * E_y) * ((K_u_ * L_u / r_g_) ** 2)


# 4.4.88
def delta(C_m_type, axial_compression_end_condition, P, F, M, L, D_o, t, E_y, S_y, S_u, material_type, M_2=0):
    C_m_ = C_m(C_m_type, M, M_2)
    f_a_ = f_a(F, D_o, D_o - 2 * t) + f_q(P, D_o, D_o - 2 * t)
    FS_ = FS_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    F_e_ = F_e(axial_compression_end_condition, D_o, t, L, E_y)

    return C_m_ / (1 - (f_a_ * FS_ / F_e_))


# 4.4.85, 86, 87
def part_h_acceptability(C_m_type, axial_compression_end_condition, P, F, M, V, phi, L, D_o, t, E_y, S_y, S_u,
                         material_type, M_2=0):
    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)
    f_a_ = f_a(F, D_o, D_o - 2 * t) + f_q(P, D_o, D_o - 2 * t)
    f_b_ = f_b(M, S_y)
    K_S_ = K_S(V, phi, L, D_o, t, E_y, S_y, S_u, material_type)
    F_xha_ = F_xha(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type)
    F_bha_ = F_bha(P, M, D_o, t, L, E_y, S_y, S_u, material_type)
    delta_ = delta(C_m_type, axial_compression_end_condition, P, F, M, L, D_o, t, E_y, S_y, S_u, material_type, M_2)

    if l_c <= 0.15:
        return (((f_a_ / (K_S_ * F_xha_)) ** 1.7) + (f_b_ / (K_S_ * F_bha_))) <= 0
    elif l_c <= 1.2:
        if (f_a_ / (K_S_ * F_xha_)) >= 0.2:
            return ((f_a_ / (K_S_ * F_xha_)) + ((8 * delta_ * f_b_) / (9 * K_S_ * F_bha_))) <= 1.0
        else:
            return ((f_a_ / (2 * K_S_ * F_xha_)) + ((delta_ * f_b_) / (K_S_ * F_bha_))) <= 1.0
    else:
        CalcExceptions.Div2CylinderException("")


# 4.4.12.2 (i), Axial Compressive Stress, Compressive Bending Stress, and Shear
def part_i_acceptability(axial_compression_end_condition, C_m_type, F, V, M, P, phi, L, D_o, t, E_y, S_y, S_u,
                         material_type, M_2=0):
    l_c = lambda_c(axial_compression_end_condition, L, D_o, D_o - 2 * t, E_y, S_y, S_u, material_type)
    f_a_ = f_a(F, D_o, D_o - 2 * t) + f_q(P, D_o, D_o - 2 * t)
    f_b_ = f_b(M, S_y)
    K_S_ = K_S(V, phi, L, D_o, t, E_y, S_y, S_u, material_type)
    F_ca_ = F_ca(axial_compression_end_condition, L, D_o, t, E_y, S_y, S_u, material_type)
    F_xa_ = F_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    F_ba_ = F_ba(D_o, t, L, E_y, S_y, S_u, material_type)
    delta_ = delta(C_m_type, axial_compression_end_condition, P, F, M, L, D_o, t, E_y, S_y, S_u, material_type, M_2)

    if l_c <= 0.15:
        return ((f_a_ / (K_S_ * F_xa_)) ** 1.7 + (f_b_ / (K_S_ * F_ba_))) <= 1.0
    elif l_c <= 1.2:
        if (f_a_ / (K_S_ * F_ca_)) >= 0.2:
            return ((f_a_ / (K_S_ * F_ca_)) + (8 * delta_ * f_b_ / (9 * K_S_ * F_ba_))) <= 1.0
        else:
            return ((f_a_ / (2 * K_S_ * F_ca_)) + (delta_ * f_b_ / (K_S_ * F_ba_))) <= 1.0
    else:
        CalcExceptions.Div2CylinderException("lambda_c is out of rang for part (i)")


# 4.4.12.2 (j), Maximum deviation modifier for 4.4.12.2


# 4.4.12.2 (k), common expressions in (a)-(j)
# 4.4.95
def A(D_o, D_i):
    return 0.25 * math.pi * (D_o ** 2 - D_i ** 2)


# 4.4.96
def S(D_o, D_i):
    return (math.pi * (D_o ** 2 - D_i ** 2)) / (32 * D_o)


# 4.4.97
def f_h(P, D_o, t):
    return (P * D_o) / (2 * t)


# 4.4.98
def f_b(M, S_):
    return M / S_


# 4.4.99
def f_a(F, D_o, D_i):
    return F / A(D_o, D_i)


# 4.4.100
def f_q(P, D_o, D_i):
    return (P * math.pi * D_i ** 2) / (4 * A(D_o, D_i))


# 4.4.101
def f_v(V, phi, D_o, D_i):
    return V * math.sin(phi) / A(D_o, D_i)


# 4.4.102
def r_g(D_o, D_i):
    return 0.25 * math.sqrt(D_o ** 2 + D_i ** 2)


# 4.4.104
def lambda_c(axial_compression_end_condition, L_u, D_o, D_i, E_y, S_y, S_u, material_type):
    F_xa_ = F_xa(D_o, (D_o - D_i) / 2, L_u, E_y, S_y, S_u, material_type)
    FS_xa_ = FS_xa(D_o, (D_o - D_i) / 2, L_u, E_y, S_y, S_u, material_type)
    K_u_ = K_u(axial_compression_end_condition)
    r_g_ = r_g(D_o, D_i)
    return (K_u_ * L_u / (math.pi * r_g_)) * math.sqrt(F_xa_ * FS_xa_ / E_y)


def get_result(axial_compression_end_condition, C_m_type, P, M, F, V, phi, D_o, t, L, E_y, S_y, S_u, material_type,
               M_2=0):
    D_i = D_o - 2 * t
    F_he_ = F_he(E_y, D_o, t, L)
    F_ic_1 = Div2Part4_4_general.F_ic(F_he_, E_y, S_y, S_u, material_type)
    S_ = S(D_o, D_i)
    F_xe_ = F_xe(E_y, L, D_o, t)
    F_ic_2 = Div2Part4_4_general.F_ic(F_xe_, E_y, S_y, S_u, material_type)
    F_xa_ = F_xa(D_o, t, L, E_y, S_y, S_u, material_type)
    F_ve_ = F_ve(L, D_o, t, E_y)
    F_ic_3 = Div2Part4_4_general.F_ic(F_ve_, E_y, S_y, S_u, material_type)

    result = {
        "M_x": M_x(L, D_o / 2, t),
        "C_h": C_h(D_o, t, L),
        "F_he": F_he_,
        "F_ic_he": F_ic_1,
        "FS_he": Div2Part4_4_general.FS(F_ic_1, S_y),
        "F_ha": F_ha(D_o, t, L, E_y, S_y, S_u, material_type),
        "P_a_": P_a(E_y, D_o, t, L, S_y, S_u, material_type),
        "A": A(D_o, D_i),
        "S": S_,
        "f_h": f_h(P, D_o, t),
        "f_b": f_b(M, S_),
        "f_a": f_a(F, D_o, D_i),
        "f_q": f_q(P, D_o, D_i),
        "f_v": f_v(V, phi, D_o, D_i),
        "r_g": r_g(D_o, D_i),
        "lambda_c": lambda_c(axial_compression_end_condition, L, D_o, D_i, E_y, S_y, S_u, material_type),
        "c_bar": c_bar(L, D_o, t),
        "C_x": C_x(L, D_o, t),
        "F_xe": F_xe_,
        "F_ic_xe": F_ic_2,
        "FS_xe": Div2Part4_4_general.FS(F_ic_2, S_y),
        "F_xa": F_xa_,
        "F_ca": F_ca(axial_compression_end_condition, L, D_o, t, E_y, S_y, S_u, material_type),
        "F_ba": F_xa_,
        "alpha_v": alpha_v(D_o, t),
        "C_v": C_v(L, D_o, t),
        "F_ve": F_ve_,
        "F_ic_ve": F_ic_3,
        "FS_ve": Div2Part4_4_general.FS(F_ic_3, S_y),
        "F_va": F_va(L, D_o, t, E_y, S_y, S_u, material_type),
        "f_x": f_x(F, P, D_o, D_i),
        "C_2": C_2(F, P, D_o, D_i),
        "C_1": C_1(D_o, t, L, E_y, S_y, S_u, material_type),
        "F_ah1": F_xha_for_e_2(axial_compression_end_condition, P, F, D_o, t, L, E_y, S_y, S_u, material_type),
        "F_ah2": F_ah2(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type),
        "F_xha": F_xha(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type),
        "F_hxa": F_hxa(axial_compression_end_condition, F, P, D_o, t, L, E_y, S_y, S_u, material_type),
        "n": n(D_o, t, L, E_y, S_y, S_u, material_type),
        "C_4": C_4(P, M, D_o, t, L, E_y, S_y, S_u, material_type),
        "C_3": C_3(P, M, D_o, t, L, E_y, S_y, S_u, material_type),
        "F_bha": F_bha(P, M, D_o, t, L, E_y, S_y, S_u, material_type),
        "F_hba": F_hba(P, M, D_o, t, L, E_y, S_y, S_u, material_type),
        "C_5": C_5(V, phi, P, D_o, t),
        "F_vha": F_vha(V, phi, P, L, D_o, t, E_y, S_y, S_u, material_type),
        "F_hva": F_hva(V, phi, P, L, D_o, t, E_y, S_y, S_u, material_type),
        "K_S": K_S(V, phi, L, D_o, t, E_y, S_y, S_u, material_type),
        "F_e": F_e(axial_compression_end_condition, D_o, t, L, E_y),
        "delta": delta(C_m_type, axial_compression_end_condition, P, F, M, L, D_o, t, E_y, S_y, S_u, material_type),
        "part_h_acceptability": part_h_acceptability(C_m_type, axial_compression_end_condition, P, F, M, V, phi, L, D_o,
                                                     t, E_y, S_y, S_u, material_type, M_2),
        "part_i_acceptability": part_i_acceptability(axial_compression_end_condition, C_m_type, F, V, M, P, phi, L, D_o,
                                                     t, E_y, S_y, S_u, material_type, M_2)
    }

    return result
