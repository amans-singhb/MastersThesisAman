##### Catalyst particle balance for WGS reactor model #####

using Pkg
Pkg.activate("WGS")
include("functions_WGSParticle.jl")

### Defining naming conventions, numeration and units ###

# for species i and j have the following numeration:
# 1 - CO
# 2 - CO2
# 3 - H2
# 4 - H2O
# 5 - N2

# For commmonly used parameters and variables:
# t - time [h]
# r - radial coordinate [m]
# T - temperature [K]
# P - pressure [atm]
# C_i - concentration of species i [mol/m^3]
# R - ideal gas constant [m^3 atm / K kmol]

#### Particle balance WGS reactor model ####

### Listing some assumptions made for this part of the code ###
# 1. T_c is kept constant, as we are only interested in the evolution of C_c_i
# 2. P_c is kept constant

### PDE system ###

# Parameters from Jacobian
D_rct_val = 4.0 # [m] 
rad_cat = 0.006 # [m]
D_cat_val = rad_cat * 2 # [m] 
d_cat_val = 5904 # [kg/m^3] [param5]

θ_val = 0.55 #[-] [param3]
τ_val = 5 #[-] [param4]

F_0 = 18398000 / 60.33863907329515 # [mol/h] (division due to inconsistency with concentrations)
T_c_val = 505.15 # [K] [param1]
P_c_val = 54.269 # [atm] [param6]
R_atmm3 = 8.2057e-2 # [m3 atm/kmol K] [param2]

C_i_val = [47.53114319; 265.2948608; 260.3569031; 381.4163208; 15.96712494;]
C_c_i_init =  C_i_val * 0.75

# Mass transfer coefficient
D_i_m_bulk = D_i_m_func(C_i_val, θ_val, τ_val, T_c_val, P_c_val)
k_c_i_val = k_c_i_func(T_c_val, P_c_val, R_atmm3, C_i_val, D_i_m_bulk, D_cat_val, D_rct_val, F_0)

using ModelingToolkit

## Parameters ##
@parameters t r T_c R θ τ d_cat P_c C_i[1:5] k_c_i[1:5]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..) D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r)
# @variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..) D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r) r_1(t, r) r_2(t, r) r_3(t, r) r_4(t, r) r_5(t, r)

C_c_i = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r)]
C_c_i_rad = [C_c_1(t, rad_cat), C_c_2(t, rad_cat), C_c_3(t, rad_cat), C_c_4(t, rad_cat), C_c_5(t, rad_cat)]
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]
# r_i = [r_1, r_2, r_3, r_4, r_5]

## Equations and Differential Equations ##

# Differential of function D_i_m
Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T_c, P_c))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

using ModelingToolkit: scalarize

eqs_Dim = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T_c, P_c))...]
# eqs_ri = [scalarize(r_i .~ r_i_func(C_c_i, d_cat, θ, P_c, T_c))...]
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) for i in 1:5]
# DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) + r_i[i] for i in 1:5]

eqs = [eqs_Dim...; DE4...]
# eqs = [eqs_Dim...; eqs_ri...; DE4...]

ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2], C_c_3(0.0, r) ~ C_c_i_init[3], C_c_4(0.0, r) ~ C_c_i_init[4], C_c_5(0.0, r) ~ C_c_i_init[5]]
BCS2 = [Dr(C_c_1(t, 0.0)) ~ 0.0, Dr(C_c_2(t, 0.0)) ~ 0.0, Dr(C_c_3(t, 0.0)) ~ 0.0, Dr(C_c_4(t, 0.0)) ~ 0.0, Dr(C_c_5(t, 0.0)) ~ 0.0]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]

# bcs = [ICS_C_c_i...; BCS2...]
bcs = [ICS_C_c_i...; BCS2...; BCS3...]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain (time is in [h])
domains = [t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]
# vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5]
params_scal = [T_c => T_c_val, R => R_atmm3, θ => θ_val, τ => τ_val, d_cat => d_cat_val, P_c => P_c_val]
params_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
params_vec_k_c_i = [k_c_i[i] => k_c_i_val[i] for i in 1:5]
params = [params_scal...; params_vec_C_i...; params_vec_k_c_i...;]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)
sol = solve(prob, KenCarp47(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
# sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
sols = sol[C_c_1(t, r)]

# Plotting 
time = 0.1
index_sol = Int(time/0.001)
solution = sols[1:index_sol, :]

using Plots

plot(solution)

using DelimitedFiles
folder = "WGS_particle/results_particle_ind"
write_to_csv("C_c_1_ind.csv", sol[C_c_1(t, r)], folder)
write_to_csv("C_c_2_ind.csv", sol[C_c_2(t, r)], folder)
write_to_csv("C_c_3_ind.csv", sol[C_c_3(t, r)], folder)
write_to_csv("C_c_4_ind.csv", sol[C_c_4(t, r)], folder)
write_to_csv("C_c_5_ind.csv", sol[C_c_5(t, r)], folder)

#### Particle balance for diffusion within a sphere ####

### Listing some assumptions made for this part of the code ###
# 1. constant diffusivity
# 2. no Reaction
# 3. constant temperature and pressure
