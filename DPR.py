import numpy as np

# Coefficients
A1 = 0.31506237
A2 = -1.0467099
A3 = -0.57832729
A4 = 0.53530771
A5 = -0.61232032
A6 = -0.10488813
A7 = 0.68157001
A8 = 0.68446549

def compute_reduced_density_DPR(Tpr, Ppr):
    if Tpr <= 0 or Ppr <= 0:
        raise ValueError("Tpr and Ppr must be positive.")
    rho = 0.27 * Ppr / Tpr
    return rho

def compute_T_function_DPR(Tpr, Ppr):
    if Tpr <= 0:
        raise ValueError("Tpr must be positive.")

    T1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3))
    T2 = A4 + (A5 / Tpr)
    T3 = A5 * A6 / Tpr
    T4 = A7 / (Tpr ** 3)
    T5 = 0.27 * Ppr / Tpr

    return T1, T2, T3, T4, T5

def compute_rho_function_DPR(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    T1, T2, T3, T4, T5 = compute_T_function_DPR(Tpr, Ppr)

    rho_function = 1 + T1 * rho + T2 * (rho ** 2) + T3 * (rho ** 5) + (
                T4 * (rho ** 2) * (1 + (A8 * (rho ** 2))) * np.exp(-A8 * (rho ** 2))) - (T5 / rho)

    return rho_function

def compute_derivative_function_DPR(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    T1, T2, T3, T4, T5 = compute_T_function_DPR(Tpr, Ppr)
    df_function = T1 + 2 * T2 * rho + 5 * T3 * rho ** 4 + \
                  2 * T4 * rho * np.exp(-A8 * rho ** 2) * (
                              (1 + 2 * A8 * rho ** 2) - A8 * rho ** 2 * (1 + A8 * rho ** 2)) + T5 / rho ** 2
    return df_function

def compute_effective_reduced_density_DPR(Tpr, Ppr, tol=1e-12, runs=200):
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")

    rho = compute_reduced_density_DPR(Tpr, Ppr)

    rho_1 = rho
    i = 1

    while i <= runs:
        if abs(compute_rho_function_DPR(Tpr, Ppr, rho_1)) < tol:
            return rho_1
        else:
            rho_1 = rho - (compute_rho_function_DPR(Tpr, Ppr, rho) / compute_derivative_function_DPR(Tpr, Ppr, rho))
            rho = rho_1
            i += 1
    return rho_1

def compute_Z_DPR(sp_gravity, temperature, pressure, fluid_type):
    if sp_gravity <= 0 or temperature <= 0 or pressure <= 0:
        raise ValueError("Specific gravity, temperature, and pressure must be positive.")

    if fluid_type.upper() == 'NATURAL GAS':
        Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
        Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)

    elif fluid_type.upper() == 'CONDENSATE':
        Tc = 187 + 330 * sp_gravity - 71.5 * (sp_gravity ** 2)
        Pc = 706 - 51.7 * sp_gravity - 11.1 * (sp_gravity ** 2)

    else:
        raise ValueError("Invalid fluid type. Must be 'NATURAL GAS' or 'CONDENSATE'.")

    Ppr = round((pressure / Pc), 2)
    Tpr = round((temperature / Tc), 2)

    if Tpr > 1:
        rho_effective = compute_effective_reduced_density_DPR(Tpr, Ppr)
        if rho_effective is None:
            return None
        T1, T2, T3, T4, T5 = compute_T_function_DPR(Tpr, Ppr)
        z = round((1 + T1 * rho_effective + T2 * (rho_effective ** 2) + T3 * (rho_effective ** 5) +
                   (T4 * (rho_effective ** 2) * (1 + (A8 * (rho_effective ** 2))) * np.exp(
                       -A8 * (rho_effective ** 2)))), 4)
        return z, Pc, Tc
    else:
        raise ValueError('Pseudo-reduced temperature must be greater than 1')
