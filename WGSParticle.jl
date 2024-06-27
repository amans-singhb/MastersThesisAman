##### Catalyst particle balance for WGS reactor model #####

using Pkg
Pkg.activate("WGS")

### Defining naming conventions and numeration ###

# for species i and j have the following numeration:
# 1 - CO
# 2 - CO2
# 3 - H2
# 4 - H2O
# 5 - N2

#### Define functions for physical properties ####

using Symbolics

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
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    
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
function μ_mix_func(y, T, M_i)
    μ_i = [μ_i_func(T, i) for i in 1:5]

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

## Other functions ##

# Reaction rate of if
function r_i_func(C_i, d_cat, θ, P, T)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]

    R_joule = 8.314
    K_eq = exp((4577.8 / T) - 4.33)
    r_co_min = d_cat * (1 - θ) * (2.96e5) * exp(-47400 / (R_joule * T)) * P^2 * (y[1] * y[4] - ((y[2] * y[3]) / K_eq))

    r_i = [-r_co_min, -r_co_min, r_co_min, r_co_min, 0]

    return r_i
end

# Bed porosity
function ϵ_b_func(D_rct, D_cat)
    0.38 + 0.073 * (1 - (D_rct / D_cat - 2)^2 / (D_rct / D_cat)^2)
end

# Molar flux
function G_func(F_0, D_rct, ϵ_b)
    (4 * F_0) / (pi * D_rct^2 * ϵ_b)
end

# Mass transfer coefficient
function k_c_i_func(T_c, P_c, R, C_i, D_i_m)
    ρ = P_c / (R * T_c)
    
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    M_i = [28.01, 44.01, 2.016, 18.016, 28.014]
    M = sum(y .* M_i)

    μ = μ_mix_func(y, T_c, M_i)

    D_cat = 0.25e-3
    D_rct = 12.7e-3 # [m]
    F = 10 # [mol/h]
    
    ϵ_b = ϵ_b_func(D_rct, D_cat)
    G = G_func(F, D_rct, ϵ_b)

    k_c_i =  0.357 * abs.(((ρ * M * D_i_m) / μ)).^(2/3) * (G / (ρ * M * ϵ_b)) * abs((μ / (D_cat * G)))^0.359
    
    return k_c_i
end

#### Particle balance WGS reactor model ####

### Listing some assumptions made for this part of the code ###
# 1. T_c is kept constant, as we are only interested in the evolution of C_c_i
# 2. P_c is kept constant

### PDE system ###

# Reactor and catalyst properties
D_cat_val = 0.25e-1 # [m]
rad_cat = 0.5 * D_cat_val # [m]
L_val = 4.8e-3 # [m]

d_cat_val = 5904 # [kg/m^3]

# Constants 
θ_val = 0.55 #[-]
τ_val = 5 #[-]

# Temperature and pressure values (constant)
T_val = 450.0 # [K]
T_c_val = 300.0 # [K]
P_c_val = 1.3 # [atm]

# Parameters for inlet and initial conditions
F_0 = 10 # [mol/h]
R_atmm3 = 8.2057e-5 # [m3 atm/mol K]

# Inlet conditions
T_in = T_val # [K] 200 C
P_in = P_c_val # [atm]
y_0 = [0.2749, 0.1198, 0.2686, 0.3165, 0.0202]
V_flow_0 = (F_0 * R_atmm3 * T_in) / P_in # [m3/h]
C_i_val = y_0 * (F_0 / V_flow_0) # [mol/m3]

# Initial conditions
T_init = T_c_val # [K] 200 C
P_init = P_c_val # [atm]
y_init = [0.0, 0.0, 0.0, 0.0, 1.0]
V_flow_init = (F_0 * R_atmm3 * T_init) / P_init # [m3/h]
C_c_i_init = y_init * (F_0 / V_flow_init) # [mol/m3]

using ModelingToolkit

## Parameters ##
@parameters t r T_c R θ τ rad d_cat P_c C_i[1:5]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..) D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r) r_1(t, r) r_2(t, r) r_3(t, r) r_4(t, r) r_5(t, r) k_c_1(t) k_c_2(t) k_c_3(t) k_c_4(t) k_c_5(t)

C_c_i = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r)]
C_c_i_rad = [C_c_1(t, rad_cat), C_c_2(t, rad_cat), C_c_3(t, rad_cat), C_c_4(t, rad_cat), C_c_5(t, rad_cat)]
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]
r_i = [r_1, r_2, r_3, r_4, r_5]
k_c_i = [k_c_1, k_c_2, k_c_3, k_c_4, k_c_5]

## Equations and Differential Equations ##

# Differential of function D_i_m
Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T_c, P_c))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

using ModelingToolkit: scalarize

#eqs_Pc = [P_c ~ (C_c_i[1] + C_c_i[2] + C_c_i[3] + C_c_i[4] + C_c_i[5]) * R * T_c]
eqs_Dim = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T_c, P_c))...]
eqs_ri = [scalarize(r_i .~ r_i_func(C_c_i, d_cat, θ, P_c, T_c))...]
eqs_kci = [scalarize(k_c_i .~ k_c_i_func(T_c, P_c, R, C_i, D_i_m))...]
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / rad) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) + r_i[i] for i in 1:5]

eqs = [eqs_Dim...; eqs_ri...; eqs_kci...; DE4...]

ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2], C_c_3(0.0, r) ~ C_c_i_init[3], C_c_4(0.0, r) ~ C_c_i_init[4], C_c_5(0.0, r) ~ C_c_i_init[5]]
BCS2 = [Dr(C_c_1(t, 0.0)) ~ 0.0, Dr(C_c_2(t, 0.0)) ~ 0.0, Dr(C_c_3(t, 0.0)) ~ 0.0, Dr(C_c_4(t, 0.0)) ~ 0.0, Dr(C_c_5(t, 0.0)) ~ 0.0]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]

bcs = [ICS_C_c_i...; BCS2...; BCS3...]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain
domains = [t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5, k_c_1, k_c_2, k_c_3, k_c_4, k_c_5]
params_scal = [T_c => T_c_val, R => R_atmm3, θ => θ_val, τ => τ_val, rad => rad_cat, d_cat => d_cat_val, P_c => P_c_val]
params_vec = [C_i[i] => C_i_val[i] for i in 1:5]
params = [params_scal...; params_vec...]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = rad_cat/10
order = 2
discretization = MOLFiniteDifference([r => dr], t, order = order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)
sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
sols = sol[C_c_1(t, r)]
plot(sols[])

#### Particle balance for diffusion within a sphere ####

### Listing some assumptions made for this part of the code ###
# 1. constant diffusivity
# 2. no Reaction
# 3. constant temperature and pressure

