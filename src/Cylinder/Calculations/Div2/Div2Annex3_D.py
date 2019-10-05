import enum
import math


class Table3_D_1_Material(enum.Enum):
    FERRITIC_STEEL = enum.auto()
    STAINLESS_STEEL_AND_NICKEL_BASE_ALLOYS = enum.auto()
    DUPLEX_STAINLESS_STEEL = enum.auto()
    PRECIPITATION_HARDENABLE_NICKEL_BASE = enum.auto()
    ALUMINUM = enum.auto()
    COPPER = enum.auto()
    TITANIUM_AND_ZIRCONIUM = enum.auto()


# Div 2 table 3-D.1
def table_3_D_1_temperature_limit(material_type):
    if material_type == Table3_D_1_Material.FERRITIC_STEEL:
        return 900.0
    elif material_type == Table3_D_1_Material.STAINLESS_STEEL_AND_NICKEL_BASE_ALLOYS:
        return 900.0
    elif material_type == Table3_D_1_Material.DUPLEX_STAINLESS_STEEL:
        return 900.0
    elif material_type == Table3_D_1_Material.PRECIPITATION_HARDENABLE_NICKEL_BASE:
        return 1000.0
    elif material_type == Table3_D_1_Material.ALUMINUM:
        return 250.0
    elif material_type == Table3_D_1_Material.COPPER:
        return 150.0
    elif material_type == Table3_D_1_Material.TITANIUM_AND_ZIRCONIUM:
        return 500.0
    else:
        ValueError("The material type is not a member of the enum class Table3_D_1_Material")


# From Div 2 table 3-D.1
def m_2(material_type, R):
    const1 = 0
    const2 = 0
    if material_type == Table3_D_1_Material.FERRITIC_STEEL:
        const1 = 0.6
        const2 = 1.0
    elif material_type == Table3_D_1_Material.STAINLESS_STEEL_AND_NICKEL_BASE_ALLOYS:
        const1 = 0.75
        const2 = 1.0
    elif material_type == Table3_D_1_Material.DUPLEX_STAINLESS_STEEL:
        const1 = 0.7
        const2 = 0.95
    elif material_type == Table3_D_1_Material.PRECIPITATION_HARDENABLE_NICKEL_BASE:
        const1 = 1.90
        const2 = 0.93
    elif material_type == Table3_D_1_Material.ALUMINUM:
        const1 = 0.52
        const2 = 0.98
    elif material_type == Table3_D_1_Material.COPPER:
        const1 = 0.5
        const2 = 1.0
    elif material_type == Table3_D_1_Material.TITANIUM_AND_ZIRCONIUM:
        const1 = 0.5
        const2 = 0.98
    else:
        ValueError("The material type is not a member of the enum class Table3_D_1_Material")

    return const1 * (const2 - R)


# From Div 2 table 3-D.1
def epsilon_p(material_type):
    if material_type == Table3_D_1_Material.FERRITIC_STEEL:
        return 0.00002
    elif material_type == Table3_D_1_Material.STAINLESS_STEEL_AND_NICKEL_BASE_ALLOYS:
        return 0.00002
    elif material_type == Table3_D_1_Material.DUPLEX_STAINLESS_STEEL:
        return 0.00002
    elif material_type == Table3_D_1_Material.PRECIPITATION_HARDENABLE_NICKEL_BASE:
        return 0.000005
    elif material_type == Table3_D_1_Material.ALUMINUM:
        return 0.000005
    elif material_type == Table3_D_1_Material.COPPER:
        return 0.000005
    elif material_type == Table3_D_1_Material.TITANIUM_AND_ZIRCONIUM:
        return 0.000005
    else:
        ValueError("The material type is not a member of the enum class Table3_D_1_Material")


# 3-D.11
def epsilon_ys():
    return 0.002


# 3-D.10
def R(sigma_ys, sigma_uts):
    return sigma_ys / sigma_uts


# 3-D.12
def K(R):
    return 1.5 * (R ** 1.5) - 0.5 * (R ** 2.5) - R ** 3.5


# 3-D.9
def H(sigma_t, sigma_ys, sigma_uts):
    K_ = K(R(sigma_ys, sigma_uts))
    return (2 * (sigma_t - (sigma_ys + K_ * (sigma_uts - sigma_ys)))) / (K_ * (sigma_uts - sigma_ys))


# 3-D.8
def A_2(sigma_ys, sigma_uts, material_type):
    R_ = R(sigma_ys, sigma_uts)
    m_2_ = m_2(material_type, R_)
    return (sigma_uts * math.e ** m_2_) / (m_2_ ** m_2_)


# 3-D.7
def epsilon_2(sigma_t, sigma_ys, sigma_uts, material_type):
    A_2_ = A_2(sigma_ys, sigma_uts, material_type)
    R_ = R(sigma_ys, sigma_uts)
    m_2_ = m_2(material_type, R_)
    return (sigma_t / A_2_) ** (1 / m_2_)


# 3-D.6
def m_1(sigma_ys, sigma_uts, material_type):
    R_ = R(sigma_ys, sigma_uts)
    e_p = epsilon_p(material_type)

    return (math.log(R_) + (e_p - epsilon_ys())) / math.log(math.log(1 + e_p) / math.log(1 + epsilon_ys()))


# 3-D.5
def A_1(sigma_ys, sigma_uts, material_type):
    m_1_ = m_1(sigma_ys, sigma_uts, material_type)
    return (sigma_ys * (1 + epsilon_ys())) / (math.log(1 + epsilon_ys()) ** m_1_)


# 3-D.4
def epsilon_1(sigma_t, sigma_ys, sigma_uts, material_type):
    A_1_ = A_1(sigma_ys, sigma_uts, material_type)
    m_1_ = m_1(sigma_ys, sigma_uts, material_type)

    return (sigma_t / A_1_) ** (1 / m_1_)


# 3-D.2
def gamma_1(sigma_t, sigma_ys, sigma_uts, material_type):
    e_1 = epsilon_1(sigma_t, sigma_ys, sigma_uts, material_type)
    H_ = H(sigma_t, sigma_ys, sigma_uts)

    return 0.5 * e_1 * (1.0 - math.tanh(H_))


# 3-D.3
def gamma_2(sigma_t, sigma_ys, sigma_uts, material_type):
    e_2 = epsilon_2(sigma_t, sigma_ys, sigma_uts, material_type)
    H_ = H(sigma_t, sigma_ys, sigma_uts)

    return 0.5 * e_2 * (1.0 + math.tanh(H_))


# 3-D.1
def epsilon_t(sigma_t, E_y, sigma_ys, sigma_uts, material_type):
    g_1 = gamma_1(sigma_t, sigma_ys, sigma_uts, material_type)
    g_2 = gamma_2(sigma_t, sigma_ys, sigma_uts, material_type)
    return (sigma_t / E_y) + g_1 + g_2


# 3-D.13
def sigma_uts_t(sigma_uts, sigma_ys, material_type):
    R_ = R(sigma_ys, sigma_uts)
    m_2_ = m_2(material_type, R_)

    return sigma_uts * math.e ** m_2_


# 3-D.17
def D_1(sigma_t, sigma_ys, sigma_uts, material_type):
    A_1_ = A_1(sigma_ys, sigma_uts, material_type)
    m_1_ = m_1(sigma_ys, sigma_uts, material_type)

    return (sigma_t ** ((1 / m_1_) - 1)) / (2 * m_1_ * A_1_ ** (1 / m_1_))


# 3-D.18
def D_2(sigma_t, sigma_ys, sigma_uts, material_type):
    A_1_ = A_1(sigma_ys, sigma_uts, material_type)
    m_1_ = m_1(sigma_ys, sigma_uts, material_type)
    K_ = K(R(sigma_ys, sigma_uts))
    H_ = H(sigma_t, sigma_ys, sigma_uts)

    term1 = - 0.5 * (1 / (A_1_ ** (1 / m_1_)))
    term2 = sigma_t ** (1 / m_1_)
    term3 = 2 / (K_ * (sigma_uts - sigma_ys))
    term4 = 1 - math.tanh(H_) ** 2
    term5 = (1 / m_1_) * (sigma_t ** ((1 / m_1_) - 1)) * math.tanh(H_)

    return term1 * (term2 * term3 * term4 + term5)


# 3-D.19
def D_3(sigma_t, sigma_ys, sigma_uts, material_type):
    A_2_ = A_2(sigma_ys, sigma_uts, material_type)
    m_2_ = m_2(material_type, R(sigma_ys, sigma_uts))

    return (sigma_t ** ((1 / m_2_) - 1)) / (2 * m_2_ * A_2_ ** (1 / m_2_))


# 3-D.18
def D_4(sigma_t, sigma_ys, sigma_uts, material_type):
    A_2_ = A_2(sigma_ys, sigma_uts, material_type)
    m_2_ = m_2(material_type, R(sigma_ys, sigma_uts))
    K_ = K(R(sigma_ys, sigma_uts))
    H_ = H(sigma_t, sigma_ys, sigma_uts)

    term1 = 0.5 * (1 / (A_2_ ** (1 / m_2_)))
    term2 = sigma_t ** (1 / m_2_)
    term3 = 2 / (K_ * (sigma_uts - sigma_ys))
    term4 = 1 - math.tanh(H_) ** 2
    term5 = (1 / m_2_) * (sigma_t ** ((1 / m_2_) - 1)) * math.tanh(H_)

    return term1 * (term2 * term3 * term4 + term5)


# 3-D.16
def E_t(sigma_t, E_y, sigma_ys, sigma_uts, material_type):
    D_1_ = D_1(sigma_t, sigma_ys, sigma_uts, material_type)
    D_2_ = D_2(sigma_t, sigma_ys, sigma_uts, material_type)
    D_3_ = D_3(sigma_t, sigma_ys, sigma_uts, material_type)
    D_4_ = D_4(sigma_t, sigma_ys, sigma_uts, material_type)

    return ((1 / E_y) + D_1_ + D_2_ + D_3_ + D_4_) ** -1
