#### WGS reactor model ####

using Pkg
Pkg.activate("WGS")
include("functions_WGS.jl")

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

### Define model equations ###

using Symbolics
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

# Catalyst properties (Cu/ZnO/Al2O3)
d_cat_val = 5904 # [kg/m^3]
ρ_cat_val = d_cat_val / (82.416095 * 0.001) # [mol/m^3]
λ_cat_val = 756 # [J/ h m K] = 0.21 J/s m K
C_p_cat_val = 35.5 # placeholder value  (must be added to variables)

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
dz = L_val/10
dr = rad_cat/10
order = 2
discretization = MOLFiniteDifference([z => dz, r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

# Solving ODE
sol = solve(prob, Tsit5(), saveat=0.1)

# Saving results
using DelimitedFiles
folder = "WGS_full/results_WGS_lab"
write_to_csv("C_c_1_WGS.csv", sol[C_c_1(t, r)], folder)
write_to_csv("C_c_2_lab.csv", sol[C_c_2(t, r)], folder)
write_to_csv("C_c_3_lab.csv", sol[C_c_3(t, r)], folder)
write_to_csv("C_c_4_lab.csv", sol[C_c_4(t, r)], folder)
write_to_csv("C_c_5_lab.csv", sol[C_c_5(t, r)], folder)
