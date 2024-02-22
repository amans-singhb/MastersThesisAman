### Define functions for physical properties ###

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

# Effective diffusivity of i in j ###[TESTED]###
function D_eff_ij_func(D_ij, θ, τ)
    D_ij * θ / τ
end

# Effective diffusivity of i in the mixture ###[TESTED]###
function D_i_m_func(y, i, D_eff_ij)
    denominator = 0
    for j in eachindex(y)
        if j != i
            denominator += y[j] / D_eff_ij[j]
        end
    end

    return (1-y[i]) / denominator
end

## Heat capacity ##

# Heat capacity for i {T [K], C_p_i [cal/mol K]}
function C_p_i_func(T, A, B, C, D)
    A + B * T + C * T^2 + D / T^2
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

# Viscosity of gas phase mixture ###[TESTED]###
function μ_mix_func(y, μ, M)
    μ_mix = 0

    for i in eachindex(y)
        numerator = y[i] * μ[i]
        denominator = 0
        for j in eachindex(y)
            denominator += y[j] * sqrt(M[j] / M[i])
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

# Binary interaction parameter A_ij ###[TESTED]###
function A_ij_func(i, j, μ, M, T, T_boil, C)
    if i == j
        return 1.0
    else
        S_i = 1.5 * T_boil[i]
        S_j = 1.5 * T_boil[j]
        S_ij = C * (S_i * S_j)^0.5

        A_ij = (1/4) * (1 + ((μ[i]/μ[j]) * (M[j]/M[i])^0.75 * ((T + S_i) / (T + S_j)) )^0.5 )^2 * ((T + S_ij) / (T + S_i))
        return A_ij
    end
end

# Thermal conductivity of mixture at atmospheric pressure ###[TESTED]###
function λ_dash_func(y, λ, μ, M, T, T_boil, C)
    λ_dash = 0

    for i in eachindex(y)
        numerator = y[i] * λ[i]
        denominator = 0
        for j in eachindex(y)
            denominator += y[j] * A_ij_func(i, j, μ, M, T, T_boil, C)
        end
        λ_dash += numerator / denominator
    end

    return λ_dash 
end

# Thermal conductivity of the bulk mixture above atmospheric pressure ###[TESTED]###
function λ_func(T_cr, P_cr, Z_cr, ρ_r, M, λ_dash)
    # Define constants A, B and C based on ρ_r
    if ρ_r < 0.5
        A = 2.702e-4
        B = 0.535
        C = -1.000
    elseif ρ_r >= 0.5 && ρ_r < 2.0 # added >= to include 0.5
        A = 2.528e-4
        B = 0.670
        C = -1.069
    elseif ρ_r >= 2.0 && ρ_r < 2.8 # added >= to include 2.0
        A = 0.574e-4
        B = 1.155
        C = 2.016
    else
        print("Error: ρ_r not in range, ρ_r = ", ρ_r)
        return
    end

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
function α_func(G)
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


### Correlations ###

## Transfer coefficients ##

# Mass transfer coefficient
function k_c_i_func(ρ, M, D_i_m, μ, G, ϵ_b, D_cat)
    0.357 * ((ρ * M * D_i_m) / μ)^(2/3) * (G / (ρ * M * ϵ_b)) * (μ / (D_cat * G))^0.359
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
using ModelingToolkit
using DifferentialEquations

## Differential ##
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
Drr = Differential(r)^2

@register_symbolic F_fr_func(G, D_cat, ρ, ϵ_b, Re) # is this needed?

## Parameters ##
@parameters begin
    z
    r
    # Gas phase species balance
    α
    k_c_i
    a_v
    C_c_i_surface
    # Gas phase momentum balance
    G
    D_cat
    ρ
    ϵ_b
    Re
    # Gas phase energy balance
    C_p # implementet as constant for now, check if it should be changed and how.
    R
    h_f
    a_v
    T_c_surface
    k_c_i
    H_c_i_surface
    H_i
    # Catalyst phase species balance
    r_i
    # D_i_m # check how to solve this analytically
    # Catalyst phase energy balance
    ρ_cat
    C_p_cat
    θ
    C_p_c_i
    λ_cat
end


## Variables ##
@variables begin
    # Gas phase species balance
    C_i(z)
    T(z)
    P(z)
    # Gas phase momentum balance
    # Gas phase energy balance
    # Catalyst phase species balance
    C_c_i(z, r)
    D_i_m # check how to solve this analytically
    # Catalyst phase energy balance
    T_c(z, r)
end
 

## Equations ##
# 1. Gas phase species balance ## (check if broadcasting is needed) ##
# 2. Gas phase momentum balance
# 3. Gas phase energy balance
# 4. Catalyst phase species balance
# 5. Catalyst phase energy balance

eqs = [Dt(C_i) ~ -α * (T/P) *  Dz(C_i) - C_i * α * ((1/P) * Dz(T) - (T/P^2) * Dz(P)) + k_c_i .* a_v .* (C_c_i_surface - C_i),
    Dz(P) ~ - F_fr_func(G, D_cat, ρ, ϵ_b, Re),
    C_p * (P / (R * T)) * Dt(T) ~ (- C_p) * G * Dz(T) + h_f * a_v * (T_c_surface - T) + a_v * sum(k_c_i .* (H_c_i_surface - H_i) .* (C_c_i_surface - C_i)),
    Dt(C_c_i) ~ (((2 * D_i_m) / r) + Dr(D_i_m)) .* Dr(C_c_i) + D_i_m .* Drr(C_c_i) + r_i,
    (1 - θ) * ρ_cat * C_p_cat * Dt(T_c) + θ * sum(C_p_c_i .* C_c_i * Dt(T_c)) ~ (((2 * λ_cat) / r) * Dr(T_c) + λ_cat * Drr(T_c)) + θ * (D_i_m .* Dr(C_c_i) .* C_p_c_i * Dr(T_c))]


## Boundary conditions ##
# 1. T at reactor inlet
# 2. C_i at reactor inlet
# 3. P at reactor inlet
# 4. dT_c/dz at catalyst center
# 5. dC_c_i/dz at catalyst center
# 6. conditions at catalyst surface
# 7. conditions at catalyst surface

bcs = [T(0) ~ T_in,
    C_i(0) ~ C_i_in,
    P(0) ~ P_in,
    Dz(T_c(z, 0)) ~ 0,
    Dz(C_c_i(z, 0)) ~ 0,
    k_c_i .* (C_c_i(z, 0.5 * D_cat) - C_i(z)) ~ (- D_i_m) .* Dr(C_c_i(z, 0.5 * D_cat)),
    h_f * (T_c(z, 0.5 * D_cat) - T(z)) + sum(H_i(z) .* k_c_i .* (C_c_i(z, 0.5 * D_cat) - C_i(z))) ~ (- λ_cat) * Dr(T_c(z, 0.5 * D_cat)) - sum(H_c_i(z, 0.5 * D_cat) .* D_i_m .* Dr(C_c_i(z, 0.5 * D_cat)))]

# ### Test functions ###

# ## Test μ_mix_func ##
# y = [0.5, 0.3, 0.2]
# μ_i= [11.1e-6, 9.0e-6, 8.0e-6]
# M = [16.04, 30.07, 44.10]

# μ_mix = μ_mix_func(y, μ_i, M)

# ## Test D_eff_ij_func and D_i_m_func ##
# theta = 0.8
# tau = 1.0
# y = [0.5, 0.3, 0.2]  
# D_ij_A = [1.0e-5, 0.8e-5, 0.6e-5]  
# D_eff_ij_A = zeros(length(D_ij_A))

# for i in eachindex(D_ij_A)
#     D_eff_ij_A[i] = D_eff_ij_func(D_ij_A[i], theta, tau)
# end
# D_eff_ij_A 
# D_i_m_A = D_i_m_func(y, 1, D_eff_ij_A)

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
# T_c = 304.2
# M = 44.01
# P_c = 7.383
# V_c = 0.0940
# Z_c = P_c * V_c / (0.008314 * T_c)
# # Z_c1 = 0.274
# k_dash = 0.0220
# V = 0.22809
# ρ = V_c/V
# # ρ1 = 0.4121

# λ_test = λ_func(T_c, P_c, Z_c1, 0.5, M, k_dash)
