#### Functions for WGS reactor model ####

using Pkg
Pkg.activate("WGS")

### Defining naming conventions and numeration ###

# for species i and j have the following numeration:
# 1 - CO
# 2 - CO2
# 3 - H2
# 4 - H2O
# 5 - N2

### Listing some assumptions made for certain elements of the code ###
# 1. For parts where there are exponents, abs() is used to avoid negative values (which cause problems due to complex numbers).
#    This assumes that the negative values should not be negative / that the negative values are not possible.

### Define functions for physical properties ###
using Symbolics

## Diffusivity ##

# Binary gas diffusivities for component pairs {T [K], P [atm], D_ij [cm^2/s]}
function D_ij_func_a(T, P, A, B, C, D, E, F)
    (A * abs(T)^B / P) * abs((log(C / abs(T))))^(-2 * D) * exp(-E / T - F / abs(T)^2)
end

function D_ij_func_b(P, B)
    B / P
end

function D_ij_func_c(T, P, A, B)
    (A * T + B) / P
end

# Matrix of all binary gas diffusivities for component pairs (SPECIFIC TO THE SYSTEM)
function D_ij_matrix_func(T, P, D_ij_matrix) 
    #p = [eq, i, j, A,      B,      C,  D,  E,   F] (C is set to 1 where it has no value, to avoid log(0) error, D accounts for lack of C in eq = "a" )
    p = ["a" 3 1 15.39e-3 1.548 0.316e8 1 -2.80 1067;
        "a" 3 2 3.14e-5 1.75 1 0 11.7 0;
        "b" 3 4 0 1.020 1 0 0 0;
        "c" 3 5 6.007e-3 -0.99311 1 0 0 0;
        "a" 1 2 3.15e-5 1.57 1 0 113.6 0;
        "a" 1 4 0.187e-5 2.072 1 0 0 0;
        "c" 1 5 0 0.322 1 0 0 0;
        "a" 2 4 9.24e-5 1.5 1 0 307.9 0;
        "a" 2 5 3.15e-5 1.57 1 0 113.6 0;
        "a" 4 5 0.187e-5 2.072 1 0 0 0;]
    
    for row in eachrow(p)
        i = row[2]
        j = row[3]

        if row[1] == "a"
            D_ij_matrix[i, j] = D_ij_func_a(T, P, row[4], row[5], row[6], row[7], row[8], row[9])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        elseif row[1] == "b"
            D_ij_matrix[i, j] = D_ij_func_b(P, row[5])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        elseif row[1] == "c"
            D_ij_matrix[i, j] = D_ij_func_c(T, P, row[4], row[5])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        else
            print("Invalid equation! Eq: ", row[1])
            return
        end
    end

    return D_ij_matrix
end

# Effective diffusivity of i in j ###[TESTED]###
function D_eff_ij_func(D_ij, θ, τ)
    D_ij .* (θ / τ)
end

# Effective diffusivity of i in the mixture ###[TESTED]###
function D_i_m_func(C_i, θ, τ, T, P)
    y = C_i ./ (C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5])

    D_ij_matrix = zeros(Num, 5, 5)
    D_i_m_vec = zeros(Num, 5)

    D_ij = D_ij_matrix_func(T, P, D_ij_matrix)
    D_eff_ij = D_eff_ij_func(D_ij, θ, τ)
    
    for i in eachindex(y)
        i = i[1]
        denominator = 0
        row = D_eff_ij[i, :]

        for j in eachindex(y)
            j = j[1]
            if j != i
                denominator += y[j] / row[j]
            end
        end

        D_i_m_vec[i] = (1-y[i]) / denominator
    end

    D_i_m_vec = D_i_m_vec * 10000 # for conversion

    return D_i_m_vec
end

## Heat capacity ##

# Heat capacity for i {T [K], C_p_i [cal/mol K]}
function C_p_i_func(T, i)
    # p = [i, A, B, C, D]
    p = [1 6.60 1.20e-3 0 0;
        2 10.34 2.74e-3 0 -1.955e5;
        3 6.62 0.81e-3 0 0;
        4 8.22 0.15e-3 1.34e-6 0;
        5 6.50 1.00e-3 0 0]

    row = p[i, :]
    A, B, C, D = row[2], row[3], row[4], row[5]

    return A + B * T + C * T^2 + D / abs(T)^2
end

# Heat capacity of mixture ###[TESTED]###
function C_p_func(y, C_p_i)
    sum(y .* C_p_i)
end

## Viscosity ##

# Viscosity of gas phase of species i {T [K], μ [N s/m^2]}
function μ_i_func(T, i)
    # p = [i, A, B, C, D]
    p = [1 1.1127e-6 0.5338 94.7 0;
        2 2.148e-6 0.46 290 0;
        3 1.797e-7 0.685 -0.59 140;
        4 1.7096e-8 1.1146 0 0;
        5 6.5592e-7 0.6081 54.714 0]

    row = p[i, :]
    A, B, C, D = row[2], row[3], row[4], row[5]

    return (A * abs(T)^B) / (1 + C / T + D / abs(T)^2)
end

# Viscosity of gas phase mixture ###[TESTED]###
function μ_mix_func(y, μ_i, M_i)
    μ_mix = 0

    for i in eachindex(y)
        i = i[1]
        numerator = y[i] * μ_i[i]
        denominator = 0
        for j in eachindex(y)
            j = j[1]
            denominator += y[j] * sqrt(M_i[j] / M_i[i])
        end
        μ_mix += numerator / denominator
    end

    μ_mix = μ_mix * 1000 # for conversion

    return μ_mix
end

## Thermal conductivity ##

# Thermal conductivity of gas phase of species i {T [K], λ [kcal/h m K]}
function λ_i_func(T, i)
    # p = [i, A, B, C, D]
    p = [1 5.1489e-4 0.6863 57.13 501.92;
        2 3.1728 -0.3838 964 1.86e6;
        3 2.2811e-3 0.7452 12 0;
        4 5.3345e-6 1.3973 0 0;
        5 2.8497e-4 0.7722 16.323 373.72]

    row = p[i, :]
    A, B, C, D = row[2], row[3], row[4], row[5]

    return (A * abs(T)^B) / (1 + C / T + D / abs(T)^2)
end

# Binary interaction parameter A_ij ###[TESTED]###
function A_ij_func(i, j, μ_i, M_i, T, T_boil, C)
    if i == j
        return 1.0
    else
        S_i = 1.5 * T_boil[i]
        S_j = 1.5 * T_boil[j]
        S_ij = C * abs((S_i * S_j))^0.5

        A_ij = (1/4) * (1 + abs(((μ_i[i]/μ_i[j]) * abs((M_i[j]/M_i[i]))^0.75 * ((T + S_i) / (T + S_j))))^0.5 )^2 * ((T + S_ij) / (T + S_i))
        return A_ij
    end
end

# Thermal conductivity of mixture at atmospheric pressure ###[TESTED]###
function λ_dash_func(y, λ, μ_i, M_i, T, T_boil, C)
    λ_dash = 0

    for i in eachindex(y)
        i = i[1]
        numerator = y[i] * λ[i]
        denominator = 0
        for j in eachindex(y)
            j = j[1]
            denominator += y[j] * A_ij_func(i, j, μ_i, M_i, T, T_boil, C)
        end
        λ_dash += numerator / denominator
    end

    return λ_dash 
end

# Thermal conductivity of the bulk mixture above atmospheric pressure ###[TESTED]###
function λ_func(y, T, P, R, M, λ_dash)
    # Define critical properties for each species in the mixture (T_cr_i [K] and P_cr_i [atm])
    T_cr_i = [133, 304.12, 32.97, 647.14, 126.21] # [K] (Taken from Wolfram Alpha)
    P_cr_i = [3.5, 7.39, 1.293, 22.406, 3.39] * 9.869 # [atm] (Taken from Wolfram Alpha, originally in MPa)

    # Use Kay's rule to calculate pseudocritical properties for the mixture (making the assumption that this will give acceptable results)
    T_cr = sum(y .* T_cr_i)
    P_cr = sum(y .* P_cr_i)

    # Calculate critical molar volume, molar volume for the mixture and the reduced density
    V_cr = R * T_cr / P_cr
    V = R * T / P
    ρ_r = V_cr / V

    # Calculate compressibility factor Z_cr for the mixture
    Z_cr = (P_cr * V_cr) / (R * T_cr)

    # # Define constants A, B and C based on ρ_r
    # if ρ_r < 0.5
    #     A = 2.702e-4
    #     B = 0.535
    #     C = -1.000
    # elseif ρ_r >= 0.5 && ρ_r < 2.0 # added >= to include 0.5
    #     A = 2.528e-4
    #     B = 0.670
    #     C = -1.069
    # elseif ρ_r >= 2.0 && ρ_r < 2.8 # added >= to include 2.0
    #     A = 0.574e-4
    #     B = 1.155
    #     C = 2.016
    # else
    #     print("Error: ρ_r not in range, ρ_r = ", ρ_r_val)
    #     return
    # end

    # just for testing purposes (need to fix boolean issue)
    A = 2.528e-4
    B = 0.670
    C = -1.069

    λ = λ_dash + (A * (exp(B * ρ_r) + C)) / (((abs(T_cr)^(1/6) * abs(M)^0.5) / abs(P_cr)^(2/3)) * Z_cr^5)
    λ = λ * 1000/3600 # for conversion

    return λ
end


### Supplementary functions for model equations ###

## Gas phase species balance ##

# Molar flux
function G_func(F_0, D_rct, ϵ_b)
    (4 * F_0) / (pi * D_rct^2 * ϵ_b)
end

# Alpha (momentum balance term)
function α_func(G, R)
    G * R
end

# Reynolds number
function Re_func(ρ, u, L, μ)
    (ρ * u * L) / μ
end

# Friction factor correlation for gas flow through packed tubular reactor
function F_fr_func(G, D_cat, ρ, ϵ_b, Re)
    (- G^2 / (ρ * D_cat)) * ((1 - ϵ_b) / ϵ_b^3) * (1.75 + 4.2 * ((1 - ϵ_b) / abs(Re)^(1/6)))   
end

# Molar density of the gas phase
function ρ_func(P, T, R)
    P / (R * T)
end

# Linear gas velocity
function u_func(α, T, P)
    α * (T / P)
end

# Molar fraction of i in the gas phase
function y_func(C_i, i)
    C_i[i] / (C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5])
end

# Heat capacity integrated for i {T [K], C_p_i [cal/mol K]}
function C_p_i_integrated_func(T, i)
    # p = [i, A, B, C, D]
    p = [1 6.60 1.20e-3 0 0;
        2 10.34 2.74e-3 0 -1.955e5;
        3 6.62 0.81e-3 0 0;
        4 8.22 0.15e-3 1.34e-6 0;
        5 6.50 1.00e-3 0 0]

    row = p[i, :]
    A, B, C, D = row[2], row[3], row[4], row[5]

    return A * T + B/2 * T^2 + C/3 * T^3 - D / T
end

# Enthalpy of i
function H_i_func(T, i)
    H_form = [-110.53, -393.51, 0, -241.83, 0] .* (1000/4.18) # [cal/mol]

    C_p_i_integrated_298 = C_p_i_integrated_func(298, i)
    C_p_i_integrated_T = C_p_i_integrated_func(T, i)

    H_i = H_form[i] + (C_p_i_integrated_T - C_p_i_integrated_298)

    return H_i
end

# Reaction rate of if
function r_i_func(y, d_cat, θ, P, T)
    R_joule = 8.314
    K_eq = exp((4577.8 / T) - 4.33)
    r_co_min = d_cat * (1 - θ) * (2.96e5) * exp(-47400 / (R_joule * T)) * P^2 * (y[1] * y[4] - ((y[2] * y[3]) / K_eq))

    r_i = [-r_co_min, -r_co_min, r_co_min, r_co_min, 0]

    return r_i
end

### Correlations ###

## Transfer coefficients ##

# Mass transfer coefficient
function k_c_i_func(ρ, M, D_i_m, μ, G, ϵ_b, D_cat)
    0.357 * abs.(((ρ * M * D_i_m) / μ)).^(2/3) * (G / (ρ * M * ϵ_b)) * abs((μ / (D_cat * G)))^0.359
end

# Heat transfer coefficient
function h_f_func(ϵ_b, C_p, G, M, μ, D_cat, λ)
    1.37 * (0.357 / ϵ_b) * ((C_p * G) / M) * abs((μ / ( D_cat * G)))^0.359 * abs(((λ * M)/ (C_p * μ)))^(2/3)
end

## Other correlations ##

# Bed porosity
function ϵ_b_func(D_rct, D_cat)
    0.38 + 0.073 * (1 - (D_rct / D_cat - 2)^2 / (D_rct / D_cat)^2)
end

# surface area
function a_v_func(ϵ_b, D_cat)
    6 * (1 - ϵ_b) / D_cat
end