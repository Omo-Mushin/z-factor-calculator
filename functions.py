import numpy as np

# Hall-Yarborough Functions
def compute_reduced_rho_HY(Ppr, Tpr):
    if Ppr <= 0 or Tpr <= 0:
        raise ValueError("Ppr and Tpr must be positive.")
    rho_1 = 0.0125 * Ppr * (1 / Tpr) * np.exp(-1.2 * (1 / Tpr - 1) ** 2)
    return rho_1

def calculate_coefficients_HY(t):
    if t <= 0:
        raise ValueError("t must be positive.")
    X1 = 0.06125 * t * np.exp(-1.2 * (1 - t) ** 2)
    X2 = t * (14.76 - 9.76 * t + 4.58 * t ** 2)
    X3 = t * (90.7 - 242.2 * t + 42.4 * t ** 2)
    X4 = 2.18 + 2.82 * t
    return X1, X2, X3, X4

def compute_reduced_rho_HY_function_HY(rho, X1, X2, X3, X4, Ppr):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")
    rho_function = (-X1 * Ppr + (rho + rho ** 2 + rho ** 3 - rho ** 4) / (1 - rho) ** 3 - X2 * rho ** 2 + X3 * rho ** X4)
    return rho_function

def compute_derivative_function_HY(rho, X2, X3, X4):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")
    df_rho_function = (1 + 4 * rho + 4 * rho ** 2 - 4 * rho ** 3 + rho ** 4) / (1 - rho) ** 4 - 2 * X2 * rho + X3 * X4 * rho ** (X4 - 1)
    return df_rho_function

def compute_effective_reduced_density_HY(Ppr, Tpr, tol=1e-13, verb=False):
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")
    t = 1 / Tpr
    X1, X2, X3, X4 = calculate_coefficients_HY(t)
    rho_1 = compute_reduced_rho_HY(Ppr, Tpr)
    delta = 1
    i = 1  # iterations

    while True:
        abs_rho_func = abs(compute_reduced_rho_HY_function_HY(rho_1, X1, X2, X3, X4, Ppr))
        if np.isnan(abs_rho_func):
            rho_1 += 0.1
        elif abs_rho_func < tol:
            break
        elif np.isnan(rho_1):
            raise ValueError("rho_1 is NA")
        else:
            rho_new = rho_1 - compute_reduced_rho_HY_function_HY(rho_1, X1, X2, X3, X4, Ppr) / compute_derivative_function_HY(rho_1, X2, X3, X4)
            delta = abs(rho_1 - rho_new)

            if delta < tol:
                break
            rho_1 = rho_new
            i += 1

    if verb:
        print(f"number of iterations: {i}\n Pseudo-reduced pressure: {Ppr}\n delta: {delta} effective reduced density: {rho_1} \n Pseudo-reduced temperature: {Tpr}")
    return rho_1, X1

def compute_Z_Hall_Yarborough(sp_gravity, temperature, pressure, fluid_type):
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
        rho_effective, X1 = compute_effective_reduced_density_HY(Ppr, Tpr)
        if rho_effective is None:
            return None
        z = round((X1 * Ppr / rho_effective), 4)
        return z, Pc, Tc
    else:
        raise ValueError('Pseudo-reduced temperature must be greater than 1')

# Dranchuk-Abou-Kaseem Functions
A1 = 0.3265
A2 = -1.0700
A3 = -0.5339
A4 = 0.01569
A5 = -0.05165
A6 = 0.5475
A7 = -0.7361
A8 = 0.1844
A9 = 0.1056
A10 = 0.6134
A11 = 0.7210

def compute_reduced_density_DAK(Tpr, Ppr):
    if Tpr <= 0 or Ppr <= 0:
        raise ValueError("Tpr and Ppr must be positive.")
    rho = 0.27 * Ppr / Tpr
    return rho

def compute_R_Values_DAK(Tpr, Ppr):
    if Tpr <= 0:
        raise ValueError("Tpr must be positive.")

    R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
    R2 = 0.27 * Ppr / Tpr
    R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
    R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
    R5 = A10 / (Tpr ** 3)

    return R1, R2, R3, R4, R5

def compute_rho_function_DAK(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    R1, R2, R3, R4, R5 = compute_R_Values_DAK(Tpr, Ppr)

    rho_function = R1 * rho - R2 / rho + R3 * rho ** 2 - R4 * rho ** 5 + \
                   R5 * rho ** 2 * (1 + A11 * rho ** 2) * np.exp(-A11 * rho ** 2) + 1

    return rho_function

def compute_derivative_function_DAK(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    R1, R2, R3, R4, R5 = compute_R_Values_DAK(Tpr, Ppr)

    derivative_rho_function = R1 + R2 / rho ** 2 + 2 * R3 * rho - 5 * R4 * rho ** 4 + \
                              2 * R5 * rho * np.exp(-A11 * rho ** 2) * \
                              ((1 + 2 * A11 * rho ** 3) - A11 * rho ** 2 * (1 + A11 * rho ** 2))

    return derivative_rho_function

def compute_effective_reduced_density_DAK(Tpr, Ppr, tol=1e-12, runs=200):
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")

    rho = compute_reduced_density_DAK(Tpr, Ppr)

    rho_1 = rho
    i = 1

    while 1:
        if abs(compute_rho_function_DAK(Tpr, Ppr, rho_1)) < tol:
            return rho_1

        rho_new = rho_1 - compute_rho_function_DAK(Tpr, Ppr, rho_1) / compute_derivative_function_DAK(Tpr, Ppr, rho_1)
        delta = abs(rho_1 - rho_new)

        if delta < tol:
            return rho_new  # Exit the loop if ideal delta found

        rho_1 = rho_new
        i += 1

        if np.isnan(rho_1) or np.isinf(rho_1):
            raise ValueError('Encountered NaN or infinity in density computation.')

    raise ValueError('Exceeded maximum iterations. No solution found.')

def compute_Z_with_DAK(sp_gravity, temperature, pressure, fluid_type):
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

    Ppr = pressure / Pc
    Tpr = temperature / Tc

    if (0.2 <= Ppr < 30) and (1 < Tpr <= 3.0):
        rho_effective = compute_effective_reduced_density_DAK(Tpr, Ppr)
        if rho_effective is None:
            return None
        Z = round((0.27 * Ppr / (rho_effective * Tpr)), 4)
        return Z, Pc, Tc
    else:
        raise ValueError('Ppr must be in the range [0.2, 30) and Tpr in the range (1, 3.0].')
