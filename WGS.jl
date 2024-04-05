### Defining nameing conventions and numeration ###

# for species i and j have the following numeration:
# 1 - CO
# 2 - CO2
# 3 - H2
# 4 - H2O
# 5 - N2


### Define functions for physical properties ###
using Symbolics
using ModelingToolkit

## Diffusivity ##

# Binary gas diffusivities for component pairs {T [K], P [atm], D_ij [cm^2/s]}
function D_ij_func_a(T, P, A, B, C, D, E, F)
    (A * T^B / P) * (log(C / T))^(-2 * D) * exp(-E / T - F / T^2)
end

function D_ij_func_b(P, B)
    B / P
end

function D_ij_func_c(T, P, A, B)
    (A * T + B) / P
end

# Matrix of all binary gas diffusivities for component pairs (SPECIFIC TO THE SYSTEM)
function D_ij_matrix_func(T, P)
    D_ij_matrix = Matrix{Num}(undef, 5, 5)

    #p = [eq, i, j, A,      B,      C,  D,  E,   F] (C is set to 1 where it has no value, to avoid log(0) error)
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
function D_i_m_func(y, D_eff_ij)
    D_i_m_vec = Array{Num}(undef, 5)
    
    for i in eachindex(y)
        denominator = 0
        index = i[1]
        row = D_eff_ij[index, :]

        for j in eachindex(y)
            if j != i
                denominator += y[j] / row[j]
            end
        end

        D_i_m_vec[i] = (1-y[i]) / denominator
    end

    return D_i_m_vec
end

## Heat capacity ##

# Heat capacity for i {T [K], C_p_i [cal/mol K]}
function C_p_i_func(T, A, B, C, D)
    A + B * T + C * T^2 + D / T^2
end

# Array of Heat capacity for all species (SPECIFIC TO THE SYSTEM)
function C_p_i_vector_func(T)
    C_p_i_vector = Array{Num}(undef, 5)
    
    # p = [i, A, B, C, D]
    p = [1 6.60 1.20e-3 0 0;
        2 10.34 2.74e-3 0 -1.955e5;
        3 6.62 0.81e-3 0 0;
        4 8.22 0.15e-3 1.34e-6 0;
        5 6.50 1.00e-3 0 0]
    
    for row in eachrow(p)
        i = row[1]
        C_p_i_vector[i] = C_p_i_func(T, row[2], row[3], row[4], row[5])
    end

    return C_p_i_vector
end

# Heat capacity of mixture ###[TESTED]###
function C_p_func(y, C_p_i)
    sum(y .* C_p_i)
end

## Viscosity ##

# Viscosity of gas phase of species i {T [K], μ [N s/m^2]}
function μ_i_func(T, A, B, C, D)
    (A * T^B) / (1 + C / T + D / T^2)
end

# Array of viscosity for all species (SPECIFIC TO THE SYSTEM)
function μ_i_vector_funct(T)
    μ_i_vector = Array{Num}(undef, 5)
    
    # p = [i, A, B, C, D]
    p = [1 1.1127e-6 0.5338 94.7 0;
        2 2.148e-6 0.46 290 0;
        3 1.797e-7 0.685 -0.59 140;
        4 1.7096e-8 1.1146 0 0;
        5 6.5592e-7 0.6081 54.714 0]
    
    for row in eachrow(p)
        i = Int(row[1])
        μ_i_vector[i] = μ_i_func(T, row[2], row[3], row[4], row[5])
    end

    return μ_i_vector
end

# Viscosity of gas phase mixture ###[TESTED]###
function μ_mix_func(y, μ_i, M_i)
    μ_mix = 0

    for i in eachindex(y)
        numerator = y[i] * μ_i[i]
        denominator = 0
        for j in eachindex(y)
            denominator += y[j] * sqrt(M_i[j] / M_i[i])
        end
        μ_mix += numerator / denominator
    end

    return μ_mix
end

## Thermal conductivity ##

# Thermal conductivity of gas phase of species i {T [K], λ [kcal/h m K]}
function λ_i_func(T, A, B, C, D)
    (A * T^B) / (1 + C / T + D / T^2)
end

# Array of thermal conductivity for all species (SPECIFIC TO THE SYSTEM)
function λ_i_vector_func(T)
    λ_i_vector = Array{Num}(undef, 5)
    
    # p = [i, A, B, C, D]
    p = [1 5.1489e-4 0.6863 57.13 501.92;
        2 3.1728 -0.3838 964 1.86e6;
        3 2.2811e-3 0.7452 12 0;
        4 5.3345e-6 1.3973 0 0;
        5 2.8497e-4 0.7722 16.323 373.72]

    for row in eachrow(p)
        i = Int(row[1])
        λ_i_vector[i] = λ_i_func(T, row[2], row[3], row[4], row[5])
    end

    return λ_i_vector
end

# Binary interaction parameter A_ij ###[TESTED]###
function A_ij_func(i, j, μ_i, M_i, T, T_boil, C)
    if i == j
        return 1.0
    else
        S_i = 1.5 * T_boil[i]
        S_j = 1.5 * T_boil[j]
        S_ij = C * (S_i * S_j)^0.5

        A_ij = (1/4) * (1 + ((μ_i[i]/μ_i[j]) * (M_i[j]/M_i[i])^0.75 * ((T + S_i) / (T + S_j)) )^0.5 )^2 * ((T + S_ij) / (T + S_i))
        return A_ij
    end
end

# Thermal conductivity of mixture at atmospheric pressure ###[TESTED]###
function λ_dash_func(y, λ, μ_i, M_i, T, T_boil, C)
    λ_dash = 0

    for i in eachindex(y)
        i = Int(i[1])
        numerator = y[i] * λ[i]
        denominator = 0
        for j in eachindex(y)
            j = Int(j[1])
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

    λ = λ_dash + (A * (exp(B * ρ_r) + C)) / (((T_cr^(1/6) * M^0.5) / P_cr^(2/3)) * Z_cr^5)

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
    (- G^2 / (ρ * D_cat)) * ((1 - ϵ_b) / ϵ_b^3) * (1.75 + 4.2 * ((1 - ϵ_b) / Re^(1/6)))   
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
function y_func(C_i)
    C_i / sum(C_i)
end

# Enthalpy of i
using Integrals, Cuba
function H_i_func(T)
    H_form = [-110.53, -393.51, 0, -241.83, 0] # [kJ/mol]

    domain = (298, T)
    prob = IntegralProblem(C_p_i_vector_func, domain)
    sol = solve(prob, CubaCuhre(); reltol=1e-6,abstol=1e-6)

    H_i = H_form + sol.u

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
    0.357 * ((ρ * M * D_i_m) / μ).^(2/3) * (G / (ρ * M * ϵ_b)) * (μ / (D_cat * G))^0.359
end

# Heat transfer coefficient
function h_f_func(ϵ_b, C_p, G, M, μ, D_cat, λ)
    1.37 * (0.357 / ϵ_b) * ((C_p * G) / M) * (μ / ( D_cat * G))^0.359 * ((λ * M)/ (C_p * μ))^(2/3)
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

## Registering symbolic functions ## (Double check if all are needed)
@register_symbolic D_ij_func_a(T, P, A, B, C, D, E, F) # added to eqs #
@register_symbolic D_ij_func_b(P, B) # added to eqs #
@register_symbolic D_ij_func_c(T, P, A, B) # added to eqs #
@register_symbolic D_ij_matrix_func(T, P) # added to eqs #
@register_symbolic D_eff_ij_func(D_ij, θ, τ) # added to eqs #
@register_symbolic D_i_m_func(y, D_eff_ij) # added to eqs #

@register_symbolic C_p_i_func(T, A, B, C, D) # added to eqs #
@register_symbolic C_p_i_vector_func(T) # added to eqs #
@register_symbolic C_p_func(y, C_p_i) # added to eqs #

@register_symbolic μ_i_func(T, A, B, C, D) # added to eqs #
@register_symbolic μ_i_vector_funct(T) # added to eqs #
@register_symbolic μ_mix_func(y, μ_i, M_i) # added to eqs #

@register_symbolic λ_i_func(T, A, B, C, D) # added to eqs #
@register_symbolic λ_i_vector_func(T) # added to eqs #
@register_symbolic A_ij_func(i, j, μ_i, M_i, T, T_boil, C) # added to eqs #
@register_symbolic λ_dash_func(y, λ_i, μ_i, M_i, T, T_boil, C) # added to eqs #
@register_symbolic λ_func(y, T, P, R, M, λ_dash) # added to eqs #

#@register_symbolic G_func(F_0, D_rct, ϵ_b) # constant
#@register_symbolic α_func(G, R) # constant
@register_symbolic Re_func(ρ, u, L, μ) # added to eqs #
@register_symbolic F_fr_func(G, D_cat, ρ, ϵ_b, Re) # added to eqs #
@register_symbolic ρ_func(P, T, R) # added to eqs #
@register_symbolic u_func(α, T, P) # added to eqs #
@register_symbolic y_func(C_i) # added to eqs #
@register_symbolic H_i_func(T) # added to eqs #
@register_symbolic r_i_func(y, d_cat, θ, P, T) # added to eqs #

@register_symbolic k_c_i_func(ρ, M, D_i_m, μ, G, ϵ_b, D_cat) # added to eqs #
@register_symbolic h_f_func(ϵ_b, C_p, G, M, μ, D_cat, λ) # added to eqs #
#@register_symbolic ϵ_b_func(D_rct, D_cat) # constant
#@register_symbolic a_v_func(ϵ_b, D_cat) # constant
 

using OrdinaryDiffEq, DomainSets, MethodOfLines
using DelimitedFiles

### Solving PDesystem ###

## Define and calculate parameters ##

# Variable
F_0 = 10 # [mol/h]

R_joule = 8.314 # [J/mol K]
R_atmL = 0.082057 # [L atm/mol K]
R_atmm3 = 8.2057e-5 # [m3 atm/mol K]
R_val = R_atmL # The R chosen to be used overall (R_joule is used in r_i_func, but is added in the function itself)

T_in = 473 # [K] 200 C
P_in = 1.3 # [atm]

y_0 = [0.2749, 0.1198, 0.2686, 0.3165, 0.0202]
V_flow_0 = (F_0 * R_atmm3 * T_in) / P_in # [m3/h]
C_i_in = y_0 * (F_0 / V_flow_0) # [mol/m3]

D_cat_val = 0.25e-3 # [m]
rad_cat_val = 0.5 * D_cat_val # [m]
D_rct = 12.7e-3 # [m]
L_val = 4.8e-3 # [m]

# Fixed
M_i_val = [28.01, 44.01, 2.016, 18.016, 28.014] # [g/mol]
θ_val = 0.55 #[-]
τ_val = 5 #[-]

T_boil_val = [81.65, 194.7, 20.35, 373, 77.36] # [K]

C_val = 1.0 # [-]

d_cat_val = 5904 # [kg/m^3]
ρ_cat_val = 2400 # plaeceholder value
C_p_cat_val = 35.5 # plaeceholder value
λ_cat_val = 0.1 # plaeceholder value

ϵ_b_val = ϵ_b_func(D_rct, D_cat)
G_val = G_func(F_0, D_rct, ϵ_b)
α_val = α_func(G, R)
a_v_val = a_v_func(ϵ_b, D_cat)

## Parameters ##
@parameters begin
    t
    z
    r

    # Gas phase species balance
    α # function, constant
    a_v # function, constant
    M_i[1:5]
    θ
    τ
    
    # Gas phase momentum balance
    G # function, constant
    D_cat
    rad_cat
    ϵ_b # function, constant
    L

    # Gas phase energy balance
    R
    T_boil[1:5]
    C
    
    # Catalyst phase species balance
    d_cat

    # Catalyst phase energy balance
    ρ_cat
    C_p_cat
    λ_cat
end

## Differential ##
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables begin
    # Gas phase species balance
    y(t, z, r)[1:5]
    C_i(..)[1:5]
    T(..)
    P(..)

    # Catalyst phase species balance
    C_c_i(..)[1:5]

    # Catalyst phase energy balance
    T_c(..)

    # Other
    M(y)
    D_ij(T, P)[1:5, 1:5]
    D_eff_ij(D_ij, θ, τ)[1:5, 1:5]
    D_i_m(y, D_eff_ij)[1:5]
    ρ(P, T, R)
    μ_i(T)[1:5]
    μ(y, μ_i, M_i)
    k_c_i(ρ, M, D_i_m, μ, G, ϵ_b, D_cat)
    u(α, T, P)
    Re(ρ, u, L, μ)
    C_p_i(T)[1:5]
    C_p(y, C_p_i)
    λ_i(T)[1:5]
    λ_dash(y, λ_i, μ, M_i, T, T_boil, C)
    λ(y, T, P, R, M, λ_dash)
    h_f(ϵ_b, C_p, G, M, μ, D_cat, λ)
    C_p_c_i(T_c(t, z, r))[1:5]
    H_i(T)[1:5]
    H_c_i_surface(T_c(t, z, r))[1:5]
    r_i(y, d_cat, θ, P, T, R)[1:5]
end

## Equations ##
# 21. Gas phase species balance ## (check if broadcasting is needed) ##
# 22. Gas phase momentum balance
# 23. Gas phase energy balance
# 24. Catalyst phase species balance
# 25. Catalyst phase energy balance

equations = [y .~ y_func(C_i(t,z)),
M ~ sum(y .* M_i),
D_ij .~ D_ij_matrix_func(T(t,z), P(z)),
D_eff_ij .~ D_eff_ij_func(D_ij, θ, τ),
D_i_m .~ D_i_m_func(y, D_eff_ij),
ρ ~ ρ_func(P(z), T(t,z), R),
μ_i .~ μ_i_vector_funct(T(t,z)),
μ ~ μ_mix_func(y, μ_i, M_i),
k_c_i .~ k_c_i_func(ρ, M, D_i_m, μ, G, ϵ_b, D_cat),
u ~ u_func(α, T(t,z), P(z)),
Re ~ Re_func(ρ, u, L, μ),
C_p_i .~ C_p_i_vector_func(T(t,z)),
C_p ~ C_p_func(y, C_p_i),
λ_i .~ λ_i_vector_func(T(t,z)),
λ_dash ~ λ_dash_func(y, λ_i, μ_i, M_i, T(t,z), T_boil, C),
λ ~ λ_func(y, T(t,z), P(z), R, M, λ_dash),
h_f ~ h_f_func(ϵ_b, C_p, G, M, μ, D_cat, λ),
C_p_c_i .~ C_p_i_vector_func(T_c(t, z, r)),
H_i .~ H_i_func(T(t,z)),
H_c_i_surface .~ H_i_func(T_c(t, z, rad_cat)),
r_i .~ r_i_func(y, d_cat, θ, P(z), T(t,z)),
Dz(P(z)) ~ - F_fr_func(G, D_cat, ρ, ϵ_b, Re),
C_p * (P(z) / (R * T(t,z))) * Dt(T(t,z)) ~ (- C_p) * G * Dz(T(t,z)) + h_f * a_v * (T_c(t, z, rad_cat) - T(t,z)) + a_v * sum(k_c_i .* (H_c_i_surface - H_i) .* (C_c_i(t, z, rad_cat) - C_i(t,z)))]

DE1 = [Dt(C_i(t,z)[i]) ~ -α * (T(t,z)/P(z)) *  Dz(C_i(t,z)[i]) - C_i(t,z)[i] * α * ((1/P(z)) * Dz(T(t,z)) - (T(t,z)/P(z)^2) * Dz(P(z))) + k_c_i * a_v * (C_c_i(t, z, rad_cat)[i] - C_i(t,z)[i]) for i in 1:5]
DE4 = [Dt(C_c_i(t,z,r)[i]) ~ (((2 * D_i_m[i]) / r) + Dr(D_i_m[i])) * Dr(C_c_i(t,z,r)[i]) + D_i_m[i] * Drr(C_c_i(t,z,r)[i]) + r_i[i] for i in 1:5]
DE5 = [(1 - θ) * ρ_cat * C_p_cat * Dt(T_c(t, z, r)) + θ * sum(C_p_c_i[i] * C_c_i(t, z, r)[i] * Dt(T_c(t, z, r))) ~ (((2 * λ_cat) / r) * Dr(T_c(t, z, r)) + λ_cat * Drr(T_c(t, z, r))) + θ * (D_i_m[i] .* Dr(C_c_i(t, z, r)[i]) .* C_p_c_i[i] * Dr(T_c(t, z, r))) for i in 1:5]

eqs = [equations; DE1; DE4; DE5]

## Boundary conditions ##
# 1. T at reactor inlet
# 2. C_i at reactor inlet
# 3. P at reactor inlet
# 4. dT_c/dz at catalyst center
# 5. dC_c_i/dz at catalyst center
# 6. conditions at catalyst surface
# 7. conditions at catalyst surface

boundaries = [T(t, 0) ~ T_in,
P(0) ~ P_in,
Dz(T_c(t, z, 0)) ~ 0]
BCS1 = [C_i(t, 0)[i] ~ C_i_in[i] for i in 1:5]
BCS2 = [Dz(C_c_i(t, z, 0)[i]) ~ 0 for i in 1:5]
BCS3 = [k_c_i * (C_c_i(t, z, rad_cat)[i] - C_i(t, z)[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i(t, z, rad_cat)[i]) for i in 1:5]
STEP1 = [H_c_i_surface[i] * D_i_m[i] * Dr(C_c_i(t, z, rad_cat)[i]) for i in 1:5] # the sum won't work without a for loop
BCS4 = [h_f * (T_c(t, z, rad_cat) - T(t, z)) + sum(H_i .* k_c_i .* (C_c_i(t, z, rad_cat) - C_i(t, z))) ~ (- λ_cat) * Dr(T_c(t, z, rad_cat)) - sum(STEP1)]

bcs = [boundaries; BCS1; BCS2; BCS3; BCS4]

# Domain
domains = [t ∈ Interval(0.0, 10.0),
z ∈ Interval(0.0, L_val),
r ∈ Interval(0.0, D_cat_val)]

# PDESystem(eqs, bcs, domains, independent_vars, dependent_vars, parameters)
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, z, r], [y[1:5], C_i(t, z)[1:5], T(t, z), P(z), C_c_i(t, z, r)[1:5], T_c(t, z, r), M, D_ij[1:5, 1:5], D_eff_ij[1:5, 1:5], D_i_m[1:5], ρ, μ_i[1:5], μ, k_c_i, u, Re, C_p_i[1:5], C_p, λ_i[1:5], λ_dash, λ, h_f, C_p_c_i[1:5], H_i[1:5], H_c_i_surface[1:5], r_i[1:5]], [α => α_val, a_v => a_v_val, M_i[1:5] => M_i_val[1:5], θ =>  θ_val, τ => τ_val, G => G_val, D_cat => D_cat_val, rad_cat => rad_cat_val, ϵ_b => ϵ_b_val, L => L_val, R => R_val, T_boil[1:5] => T_boil_val[1:5], C => C_val, d_cat => d_cat_val, ρ_cat => ρ_cat_val, C_p_cat => C_p_cat_val, λ_cat => λ_cat_val])

# Discretization
dz = L_val/100
dr = D_cat_val/100
order = 2
discretization = MOLFiniteDifference([z => dz, r => dr], t)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

# Solving ODE
sol = solve(prob, Tsit5(), saveat=0.1)

# Saving results
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


