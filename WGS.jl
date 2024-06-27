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


### Define model equations ###

using ModelingToolkit

### Solving PDesystem ###

## Define and calculate parameters ##

# Variable
F_0 = 10 # [mol/h]

R_joule = 8.314 # [J/mol K]
R_atmL = 0.082057 # [L atm/mol K]
R_atmm3 = 8.2057e-5 # [m3 atm/mol K]
R_val = R_atmm3 # The R chosen to be used overall (R_joule is used in r_i_func, but is added in the function itself)

# Inlet conditions
T_in = 473 # [K] 200 C
P_in = 1.3 # [atm]

y_0 = [0.2749, 0.1198, 0.2686, 0.3165, 0.0202]
V_flow_0 = (F_0 * R_atmm3 * T_in) / P_in # [m3/h]
C_i_in = y_0 * (F_0 / V_flow_0) # [mol/m3]

# Initial conditions
T_init = 473 # [K] 200 C
P_init = 1.3 # [atm]

y_init = [0.0, 0.0, 0.0, 0.0, 1.0]
V_flow_init = (F_0 * R_atmm3 * T_init) / P_init # [m3/h]
C_i_init = y_init * (F_0 / V_flow_init) # [mol/m3]

D_cat_val = 0.25e-3 # [m]
rad_cat = 0.5 * D_cat_val # [m]
D_rct = 12.7e-3 # [m]
rad_rct = 0.5 * D_rct # [m]
L_val = 4.8e-3 # [m]

# Fixed
M_i_val = [28.01, 44.01, 2.016, 18.016, 28.014] # [g/mol]
θ_val = 0.55 #[-]
τ_val = 5 #[-]

T_boil_val = [81.65, 194.7, 20.35, 373, 77.36] # [K]

C_val = 1.0 # [-]

ϵ_b_val = ϵ_b_func(D_rct, D_cat_val)
G_val = G_func(F_0, D_rct, ϵ_b_val)
α_val = α_func(G_val, R_val)
a_v_val = a_v_func(ϵ_b_val, D_cat_val)

# Catalyst properties
d_cat_val = 5904 # [kg/m^3]
ρ_cat_val = 2400 # plaeceholder value
C_p_cat_val = 35.5 # plaeceholder value
λ_cat_val = 0.1 # plaeceholder value

## Parameters ##
@parameters begin
    t
    z
    r

    # Gas phase species balance
    α # function, constant
    a_v # function, constant
    
    θ
    τ
    
    # Gas phase momentum balance
    G # function, constant
    D_cat
    # rad_cat
    ϵ_b # function, constant
    L

    # Gas phase energy balance
    R
    
    C
    
    # Catalyst phase species balance
    d_cat

    # Catalyst phase energy balance
    ρ_cat
    C_p_cat
    λ_cat
    radius_cat

    M_i[1:5]
    T_boil[1:5]
end

## Differential ##
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables begin
    # Gas phase species balance
    C_1(..)
    C_2(..)
    C_3(..)
    C_4(..)
    C_5(..)

    T_1(..)
    T_2(..)
    T_3(..)
    T_4(..)
    T_5(..)

    P(..)

    # Catalyst phase species balance
    C_c_1(..)
    C_c_2(..)
    C_c_3(..)
    C_c_4(..)
    C_c_5(..)

    # Catalyst phase energy balance
    T_c(..)

    P_c(..)

    # Other
    y_1(t, z)
    y_2(t, z)
    y_3(t, z)
    y_4(t, z)
    y_5(t, z)

    y_c_1(t, z, r)
    y_c_2(t, z, r)
    y_c_3(t, z, r)
    y_c_4(t, z, r)
    y_c_5(t, z, r)
    
    M(t, z)

    D_1_m(t, z, r)
    D_2_m(t, z, r)
    D_3_m(t, z, r)
    D_4_m(t, z, r)
    D_5_m(t, z, r)

    ρ(t, z)
    
    μ_1(t, z)
    μ_2(t, z)
    μ_3(t, z)
    μ_4(t, z)
    μ_5(t, z)

    μ(t, z)

    k_c_1(t, z)
    k_c_2(t, z)
    k_c_3(t, z)
    k_c_4(t, z)
    k_c_5(t, z)

    u(t,z)

    Re(t, z)

    C_p_1(t, z)
    C_p_2(t, z)
    C_p_3(t, z)
    C_p_4(t, z)
    C_p_5(t, z)

    C_p(t, z)
    
    λ_1(t, z)
    λ_2(t, z)
    λ_3(t, z)
    λ_4(t, z)
    λ_5(t, z)

    λ_dash(t, z)

    λ(t, z)

    h_f(t, z)
    
    C_p_c_1(t, z, r)
    C_p_c_2(t, z, r)
    C_p_c_3(t, z, r)
    C_p_c_4(t, z, r)
    C_p_c_5(t, z, r)
    
    H_1(t, z)
    H_2(t, z)
    H_3(t, z)
    H_4(t, z)
    H_5(t, z)
    
    H_c_1_surface(t, z, r)
    H_c_2_surface(t, z, r)
    H_c_3_surface(t, z, r)
    H_c_4_surface(t, z, r)
    H_c_5_surface(t, z, r)
    
    r_1(t, z, r)
    r_2(t, z, r)
    r_3(t, z, r)
    r_4(t, z, r)
    r_5(t, z, r)
end


C_i = [C_1(t, z), C_2(t, z), C_3(t, z), C_4(t, z), C_5(t, z)]
T = [T_1(t, z), T_2(t, z), T_3(t, z), T_4(t, z), T_5(t, z)]
C_c_i = [C_c_1(t, z, r), C_c_2(t, z, r), C_c_3(t, z, r), C_c_4(t, z, r), C_c_5(t, z, r)]
C_c_i_rad = [C_c_1(t, z, rad_cat), C_c_2(t, z, rad_cat), C_c_3(t, z, rad_cat), C_c_4(t, z, rad_cat), C_c_5(t, z, rad_cat)]
y = [y_1, y_2, y_3, y_4, y_5]
y_c = [y_c_1, y_c_2, y_c_3, y_c_4, y_c_5] 
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m] 
μ_i = [μ_1, μ_2, μ_3, μ_4, μ_5]
k_c_i = [k_c_1, k_c_2, k_c_3, k_c_4, k_c_5]
C_p_i = [C_p_1, C_p_2, C_p_3, C_p_4, C_p_5]
λ_i = [λ_1, λ_2, λ_3, λ_4, λ_5]
C_p_c_i = [C_p_c_1, C_p_c_2, C_p_c_3, C_p_c_4, C_p_c_5]
H_i = [H_1, H_2, H_3, H_4, H_5]
H_c_i_surface = [H_c_1_surface, H_c_2_surface, H_c_3_surface, H_c_4_surface, H_c_5_surface]
r_i = [r_1, r_2, r_3, r_4, r_5]


## Equations ##
# DE1. Gas phase species balance ## (check if broadcasting is needed) ##
# DE2. Gas phase momentum balance
# DE3. Gas phase energy balance
# DE4. Catalyst phase species balance
# DE5. Catalyst phase energy balance

import ModelingToolkit: scalarize

equations_y = [y[i] ~ y_func(C_i, i) for i in 1:5]
equations_y_c = [y_c[i] ~ y_func(C_c_i, i) for i in 1:5]
equations_μ = [μ_i[i] ~ μ_i_func( T[1], i) for i in 1:5]
equations_C_p_i = [C_p_i[i] ~ C_p_i_func( T[1], i) for i in 1:5]
equations_λ = [λ_i[i] ~ λ_i_func( T[1], i) for i in 1:5]
equations_C_p_c_i = [C_p_c_i[i] ~ C_p_i_func(T_c(t, z, r), i) for i in 1:5]
equations_H_i = [H_i[i] ~ H_i_func( T[1], i) for i in 1:5]
equations_H_c_i_surface = [H_c_i_surface[i] ~ H_i_func(T_c(t, z, rad_cat), i) for i in 1:5]
equations_P_c = [P_c(t, z, r) ~ (C_c_i[1] + C_c_i[2] + C_c_i[3] + C_c_i[4] + C_c_i[5]) * R_val * T_c(t, z, r)]
equations1 = [equations_y; equations_y_c; equations_μ; equations_C_p_i; equations_λ; equations_C_p_c_i; equations_H_i; equations_H_c_i_surface; equations_P_c]


# equations_D_ij = [scalarize(D_ij[1:5,1:5] .~ D_ij_matrix_func(T(t, z), P(t, z), matrix_D_ij))...]
# equations_D_eff_ij = [scalarize(D_eff_ij .~ D_eff_ij_func(D_ij, θ, τ))...]
equations_D_i_m = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T_c(t, z, r), P_c(t, z, r)))...]
equations_k_c_i = [scalarize(k_c_i .~ k_c_i_func(ρ, M, D_i_m, μ, G, ϵ_b, D_cat))...]
equations_r_i = [scalarize(r_i .~ r_i_func(y_c, d_cat, θ, P_c(t, z, r), T_c(t, z, r)))...]
equations2 = [equations_D_i_m; equations_k_c_i; equations_r_i]

equations3 = [M ~ sum(y .* M_i),
    ρ ~ ρ_func(P(t, z),  T[1], R),
    μ ~ μ_mix_func(y, μ_i, M_i),
    u ~ u_func(α,  T[1], P(t, z)),
    Re ~ Re_func(ρ, u, L, μ),
    C_p ~ C_p_func(y, C_p_i),
    λ_dash ~ λ_dash_func(y, λ_i, μ_i, M_i,  T[1], T_boil, C),
    λ ~ λ_func(y,  T[1], P(t, z), R, M, λ_dash),
    h_f ~ h_f_func(ϵ_b, C_p, G, M, μ, D_cat, λ)]

equations = [equations1; equations2; equations3]

Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T_c(t, z, r), P_c(t, z, r)))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

DE1 = [Dt(C_i[i]) ~ -α * ( T[i]/P(t, z)) *  Dz(C_i[i]) - C_i[i] * α * ((1/P(t, z)) * Dz( T[i]) - ( T[i]/P(t, z)^2) * Dz(P(t, z))) + k_c_i[i] * a_v * (C_c_i_rad[i] - C_i[i]) for i in 1:5]
DE2 = Dz(P(t, z)) ~ - F_fr_func(G, D_cat, ρ, ϵ_b, Re)
DE3 = C_p * (P(t, z) / (R *  T[1])) * Dt( T[1]) ~ (- C_p) * G * Dz( T[1]) + h_f * a_v * (T_c(t, z, rad_cat) -  T[1]) + a_v * sum(k_c_i .* (H_c_i_surface - H_i) .* (C_c_i_rad - C_i))
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) + r_i[i] for i in 1:5]
STEP_DE = [C_p_c_i[i] * C_c_i[i] * Dt(T_c(t, z, r)) for i in 1:5]
DE5 = [(1 - θ) * ρ_cat * C_p_cat * Dt(T_c(t, z, r)) + θ * sum(STEP_DE) ~ (((2 * λ_cat) / radius_cat) * Dr(T_c(t, z, r)) + λ_cat * Drr(T_c(t, z, r))) + θ * (D_i_m[i] * Dr(C_c_i[i]) * C_p_c_i[i] * Dr(T_c(t, z, r))) for i in 1:5]

diffequations = [DE1; DE2; DE3; DE4; DE5]

eqs = [equations; diffequations]

## Initial conditions ##
# ICS_T. T at t = 0
# ICS_P. P at t = 0
# ICS_C_i. C_i at t = 0

C_i_init_sym = [C_1(0.0, z), C_2(0.0, z), C_3(0.0, z), C_4(0.0, z), C_5(0.0, z)]
T_init_sym = [T_1(0.0, z), T_2(0.0, z), T_3(0.0, z), T_4(0.0, z), T_5(0.0, z)]

ICS_T = [T_init_sym[i] ~ T_init for i in 1:5]
ICS_P = [P(0, z) ~ P_init]
ICS_C_i = [C_i_init_sym[i] ~ C_i_init[i] for i in 1:5]
ics = [ICS_T; ICS_P; ICS_C_i]

## Boundary conditions ##
# BCS_T. T at reactor inlet
# BCS_P. P at reactor inlet
# BCS_Tc. dT_c/dz at catalyst center
# BCS1. C_i at reactor inlet
# BCS2. dC_c_i/dz at catalyst center
# BCS3. conditions at catalyst surface
# BCS4. conditions at catalyst surface

C_i_bound_sym = [C_1(t, 0.0), C_2(t, 0.0), C_3(t, 0.0), C_4(t, 0.0), C_5(t, 0.0)]
C_c_i_bound_sym = [C_c_1(t, z, 0.0), C_c_2(t, z, 0.0), C_c_3(t, z, 0.0), C_c_4(t, z, 0.0), C_c_5(t, z, 0.0)]
T_bound_sym = [T_1(t, 0.0), T_2(t, 0.0), T_3(t, 0.0), T_4(t, 0.0), T_5(t, 0.0)]

BCS_T = [T_bound_sym[i] ~ T_in for i in 1:5]
BCS_P = [P(t, 0.0) ~ P_in]
BCS_Tc = [Dr(T_c(t, z, 0)) ~ 0]
BCS1 = [C_i_bound_sym[i] ~ C_i_in[i] for i in 1:5]
BCS2 = [Dr(C_c_i_bound_sym[i]) ~ 0.0 for i in 1:5]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]
STEP_BCS = [H_c_i_surface[i] * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5] # the sum won't work without a for loop
BCS4 = [h_f * (T_c(t, z, rad_cat) -  T[i]) + sum(H_i .* k_c_i .* (C_c_i_rad - C_i)) ~ (- λ_cat) * Dr(T_c(t, z, rad_cat)) - sum(STEP_BCS) for i in 1:5]

bcs = [ics; BCS_T; BCS_P; BCS_Tc; BCS1; BCS2; BCS3; BCS4]

using OrdinaryDiffEq, DomainSets, MethodOfLines
# Domain
domains = [t ∈ Interval(0.0, 1.0),
    z ∈ Interval(0.0, L_val),
    r ∈ Interval(0.0, rad_cat)]    

vars_imp = reduce(vcat, [[C_i[i] for i in 1:5],  [T[i] for i in 1:5], P(t, z), [C_c_i[i] for i in 1:5], T_c(t, z, r), P_c(t, z, r)])
vars_other = reduce(vcat, [[y[i] for i in 1:5], [y_c[i] for i in 1:5], M, [D_i_m[i] for i in 1:5], ρ, [μ_i[i] for i in 1:5], μ, [k_c_i[i] for i in 1:5], u, Re, [C_p_i[i] for i in 1:5], C_p, [λ_i[i] for i in 1:5], λ_dash, λ, h_f, [C_p_c_i[i] for i in 1:5], [H_i[i] for i in 1:5], [H_c_i_surface[i] for i in 1:5], [r_i[i] for i in 1:5]])
vars = [vars_imp; vars_other]

params_vec_M_i = [M_i[i] => M_i_val[i] for i in 1:5]
params_vec_T_boil = [T_boil[i] => T_boil_val[i] for i in 1:5]
params_scal = [α => α_val, a_v => a_v_val,  θ =>  θ_val, τ => τ_val, G => G_val, D_cat => D_cat_val, ϵ_b => ϵ_b_val, L => L_val, R => R_val, C => C_val, d_cat => d_cat_val, ρ_cat => ρ_cat_val, C_p_cat => C_p_cat_val, λ_cat => λ_cat_val, radius_cat => rad_cat]
params = [params_scal; params_vec_M_i; params_vec_T_boil]

# PDESystem(eqs, bcs, domains, independent_vars, dependent_vars, parameters)
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, z, r], vars, params)

# Discretization
dz = round(L_val/10, sigdigits=4)
dr = round(rad_cat/10, sigdigits=4)
# dz = L_val/10
# dr = rad_cat/10
order = 2
discretization = MOLFiniteDifference([z => dz, r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

# Solving ODE
sol = solve(prob, Tsit5(), saveat=0.1)

# Saving results
using DelimitedFiles
writedlm("WGS_results.csv", sol, ",")



# ### Test functions ###

# ## Test μ_mix_func ##
# y = [0.5, 0.3, 0.2]
# μ_i= [11.1e-6, 9.0e-6, 8.0e-6]
# M = [16.04, 30.07, 44.10]

# μ_mix = μ_mix_func(y, μ_i, M)

# ## Test D_eff_ij_func and D_i_m_func ##
# theta = 0.8
# tau = 1.0
# y = [0.2, 0.2, 0.2, 0.2, 0.2]
# D_ij_mat = D_ij_matrix_func(300, 1)

# #D_eff_ij = zeros(size(D_ij_mat))

# D_eff_ij = D_eff_ij_func(D_ij_mat, theta, tau)

# D_i_m = D_i_m_func(y, D_eff_ij)


# ## Test C_p_func ##
# y = [0.5, 0.3, 0.2]
# C_p_i = [28.1, 33.2, 42.0] 

# C_p_mix = C_p_func(y, C_p_i)

# ## Test A_ij_func, λ_dash_func ##
# y_list = [0.23, 0.77]
# T_example = 373.15
# k_list = [0.02504, 0.01587]
# μ_list = [1.161e-5, 1.361e-5]
# M_list = [46.07, 50.49]
# T_boil_list = [248.3, 248.9]

# A_11 = A_ij_func(1, 1, μ_list, M_list, T_example, T_boil_list, 1.0)
# A_12 = A_ij_func(1, 2, μ_list, M_list, T_example, T_boil_list, 1.0)
# A_21 = A_ij_func(2, 1, μ_list, M_list, T_example, T_boil_list, 1.0)
# A_22 = A_ij_func(2, 2, μ_list, M_list, T_example, T_boil_list, 1.0)

# λ_dash = λ_dash_func(y_list, k_list, μ_list, M_list, T_example, T_boil_list, 1.0)

# ## Test λ_func ## 
# M = 44.01
# P_c = 7.383
# V_c = 0.0940
# Z_c = P_c * V_c / (0.008314 * T_c)
# # Z_c1 = 0.274
# k_dash = 0.0220
# V = 0.22809
# ρ = V_c/V
# # ρ1 = 0.4121

# λ_test = λ_func(y, T, P, R, M, k_dash)

# ## Test syntax/use of symbolic functions with array as ouput ##
# using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets

# function test(u)
#     return u ./ 5
# end


# @register_symbolic test(u)

# n_comp = 2
# # @register_array_symbolic test(u) begin
# #     size=(n_comp,)
# #     eltype=eltype(Vector{Num})
# # end

# # Parameters, variables, and derivatives

# @parameters t, x, p[1:n_comp], q[1:n_comp] , f
# @variables T[1:n_comp], u(..)[1:n_comp]
# Dt = Differential(t)
# Dx = Differential(x)
# Dxx = Differential(x)^2
# params = Symbolics.scalarize(reduce(vcat,[p .=> [1.5, 2.0], q .=> [1.2, 1.8], f => 2]))
# # 1D PDE and boundary conditions

# #eqs = [Dt(u(t, x)[i]) ~ p[i] * Dxx(u(t, x)[i]) for i in 1:n_comp]
# #eqs1 = [T[i] ~ test(u(t, x))[i] for i in 1:n_comp]
# eqs = [Dt(u(t, x)[i]) ~ f * p[i] * test(u(t, x)[i]) * Dxx(u(t, x)[i]) for i in 1:n_comp]
# #eqs = [eqs1; eqs2]
# bcs = [[u(0, x)[i] ~ q[i] * cos(x),
#         u(t, 0)[i] ~ sin(t),
#         u(t, 1)[i] ~ exp(-t) * cos(1),
#         Dx(u(t,0)[i]) ~ 0.0] for i in 1:n_comp]
# bcs_collected = reduce(vcat, bcs)

# # Space and time domains
# domains = [t ∈ Interval(0.0, 1.0),
#            x ∈ Interval(0.0, 1.0)]

# # PDE system

# @named pdesys = PDESystem(eqs, bcs_collected, domains, [t, x], [u(t, x)[1:n_comp]], params)


# # Method of lines discretization
# dx = 0.1
# order = 2
# discretization = MOLFiniteDifference([x => dx], t; approx_order = order)

# # Convert the PDE problem into an ODE problem
# prob = discretize(pdesys,discretization) #error occurs here

# # Solve ODE problem
# sol = solve(prob, Tsit5(), saveat=0.2)

## for derivatives of functions:
#= @parameters r, t, z
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
@variables P(t, r, z), T(z, t, r), C(r, z), (C_i(z, t))[1:5]

teste = Dr.(D_ij_matrix_func(T, P, matrix_D_ij))
derivatives = expand_derivatives.(teste)
### Register symbolic functions ###
derivatives[1, 2] =#