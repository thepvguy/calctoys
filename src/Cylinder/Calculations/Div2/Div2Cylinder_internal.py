import enum
import math

from . import CalcExceptions


# Internal pressure


# Nomenclature definitions
class ShellType(enum.Enum):
    CYLINDER = enum.auto()
    SPHERICAL = enum.auto()
    CONICAL = enum.auto()


# Hoop stress
# 4.3.3.1 Required Thickness
def t_min(P, S_, E, D_i):
    return 0.5 * D_i * (math.exp(P / (S_ * E)) - 1)


# Combined loads, 4.3.10
# 4.3.32, 35, 38
def sigma_theta_m(shell_type, P, E, D, D_o, alpha=0):
    if shell_type == ShellType.CYLINDER:
        return (P * D) / (E * (D_o - D))
    elif shell_type == ShellType.SPHERICAL:
        return (P * D ** 2) / (E * (D_o ** 2 - D ** 2))
    elif shell_type == ShellType.CONICAL:
        return (P * D) / (E * (D_o - D) * math.cos(alpha))
    else:
        CalcExceptions.Div2CylinderException("Invalid shell type")


# 4.3.33, 36, 39
def sigma_sm(shell_type, M, F, P, E, D, D_o, theta, phi=0, alpha=0):
    if shell_type == ShellType.CYLINDER:
        return (1 / E) * (((P * D ** 2) / (D_o ** 2 - D ** 2)) + (4 * F / (math.pi * (D_o ** 2 - D ** 2))) + (
                    (32 * M * D_o * math.cos(theta)) / (math.pi * (D_o ** 4 - D ** 4))))
    elif shell_type == ShellType.SPHERICAL:
        return (1 / E) * (((P * D ** 2) / (D_o ** 2 - D ** 2)) + (
                    (4 * F) / (math.pi * (D_o ** 2 - D ** 2) * math.sin(phi) ** 2)) + (
                                      (32 * M * D_o * math.cos(phi)) / (
                                          math.pi * (D_o ** 4 - D ** 4) * math.sin(phi) ** 3)))
    elif shell_type == ShellType.CONICAL:
        return (1 / E) * (((P * D ** 2) / ((D_o ** 2 - D ** 2) * math.cos(alpha))) + (
                    (4 * F) / (math.pi * (D_o ** 2 - D ** 2) * math.cos(alpha))) + ((32 * M * D_o * math.cos(theta)) / (
                    math.pi * (D_o ** 4 - D ** 4) * math.cos(alpha))))
    else:
        CalcExceptions.Div2CylinderException("Invalid shell type")


# 4.3.34, 37, 40
def tau(shell_type, M, M_t, D_o, D, theta, phi=0, alpha=0):
    if shell_type == ShellType.CYLINDER:
        return (16 * M_t * D_o) / (math.pi * (D_o ** 4 - D ** 4))
    elif shell_type == ShellType.SPHERICAL:
        return (((32 * M * D_o) / (math.pi * (D_o ** 4 - D ** 4) * math.sin(phi) ** 3)) * (
                    (math.cos(phi)) / (math.sin(phi) ** 3))) * math.sin(theta) + (
                           (16 * M_t * D_o) / (math.pi * (D_o ** 4 - D ** 4) * math.sin(phi) ** 2))
    elif shell_type == ShellType.CONICAL:
        return (32 * M * D_o / (math.pi * (D_o ** 4 - D ** 4))) * (math.tan(alpha) * math.sin(theta)) + (
                    16 * M_t * D_o / (math.pi * (D_o ** 4 - D ** 4)))
    else:
        CalcExceptions.Div2CylinderException("Invalid shell type")


# 4.3.41
def sigma_1(shell_type, M, M_t, F, P, E, D, D_o, theta, phi=0, alpha=0):
    sigma_theta_m_ = sigma_theta_m(shell_type, P, E, D, D_o, alpha)
    sigma_sm_ = sigma_sm(shell_type, M, F, P, E, D, D_o, theta, phi, alpha)
    tau_ = tau(shell_type, M, M_t, D_o, D, theta, phi, alpha)

    return 0.5 * (sigma_theta_m_ + sigma_sm_ + math.sqrt((sigma_theta_m_ - sigma_sm_) ** 2 - 4 * tau_ ** 2))


# 4.3.42
def sigma_2(shell_type, M, M_t, F, P, E, D, D_o, theta, phi=0, alpha=0):
    sigma_theta_m_ = sigma_theta_m(shell_type, P, E, D, D_o, alpha)
    sigma_sm_ = sigma_sm(shell_type, M, F, P, E, D, D_o, theta, phi, alpha)
    tau_ = tau(shell_type, M, M_t, D_o, D, theta, phi, alpha)

    return 0.5 * (sigma_theta_m_ + sigma_sm_ - math.sqrt((sigma_theta_m_ - sigma_sm_) ** 2 - 4 * tau_ ** 2))


# 4.3.43
def sigma_3():
    return 0


def sigma_vm(shell_type, M, M_t, F, P, E, D, D_o, theta, phi=0, alpha=0):
    s1 = sigma_1(shell_type, M, M_t, F, P, E, D, D_o, theta, phi, alpha)
    s2 = sigma_2(shell_type, M, M_t, F, P, E, D, D_o, theta, phi, alpha)
    s3 = sigma_3()

    return (1 / math.sqrt(2)) * math.sqrt((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2)


def compressive_stress_ok(shell_type, M, M_t, F, P, E, D, D_o, theta, S, phi=0, alpha=0):
    return sigma_vm(shell_type, M, M_t, F, P, E, D, D_o, theta, phi, alpha) <= S
