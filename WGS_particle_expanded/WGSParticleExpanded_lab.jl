##### Catalyst particle balance for WGS reactor model #####

using Pkg
Pkg.activate("WGS")
include("functions_WGSParticleExpanded.jl")

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
# 1. T is kept constant, as we are only interested in the evolution of C_c_i and T_c
# 2. P_c = P, and is kept constant

### PDE system ###

## Parameters

# Dimensions
D_rct_val = 12.7e-3 # [m] 
rad_cat = 0.125e-3 # [m]
D_cat_val = rad_cat * 2 # [m]

# Catalyst properties
d_cat_val = 5904 # [kg/m^3] [param6]
ρ_cat_val = d_cat_val / (82.416095 * 0.001) # [mol/m^3] [param7]
λ_cat_val = 756 # [J/ h m K] = 0.21 J/s m K [param8]

# Coefficients
θ_val = 0.55 #[-] [param4]
τ_val = 5 #[-] [param5]
R_atmm3 = 8.2057e-2 # [m3 atm/kmol K] [param3]

# Inlet values
F_0 = 1e-5 # [mol/h] 
T_val = 500 # [K] [param1]
P_val = 1.3 # [atm] [param2]

C_i_val = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]
C_c_i_init = C_i_val * 0.75

# Mass transfer coefficients (bulk phase)
D_i_m_bulk = D_i_m_func(C_i_val, θ_val, τ_val, T_val, P_val) # [m^2/h]
k_c_i_val = k_c_i_func(T_val, P_val, R_atmm3, C_i_val, D_i_m_bulk, D_cat_val, D_rct_val, F_0)
h_f_val = h_f_func(T_val, C_i_val, D_cat_val, D_rct_val, F_0) # [J/h m^2 K] [param9]

# Enthalpy
H_i_val = H_i_func(T_val)

using ModelingToolkit

## Parameters ##
@parameters t r T P R θ τ d_cat ρ_cat λ_cat hf
@parameters C_i[1:5] k_c_i[1:5] H_i[1:5]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..) T_c(..)
@variables D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r) r_1(t, r) r_2(t, r) r_3(t, r) r_4(t, r) r_5(t, r)
@variables C_p_c_1(t, r) C_p_c_2(t, r) C_p_c_3(t, r) C_p_c_4(t, r) C_p_c_5(t, r) H_c_1(t, r) H_c_2(t, r) H_c_3(t, r) H_c_4(t, r) H_c_5(t, r) C_p_cat(t, r)

C_c_i = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r)]
C_c_i_rad = [C_c_1(t, rad_cat), C_c_2(t, rad_cat), C_c_3(t, rad_cat), C_c_4(t, rad_cat), C_c_5(t, rad_cat)]
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]
r_i = [r_1, r_2, r_3, r_4, r_5]

C_p_c_i = [C_p_c_1, C_p_c_2, C_p_c_3, C_p_c_4, C_p_c_5]
H_c_i = [H_c_1, H_c_2, H_c_3, H_c_4, H_c_5]

## Equations and Differential Equations ##

# Differential of function D_i_m
Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T_c(t, r), P))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

using ModelingToolkit: scalarize

# Intermediate variables
Hidr = [scalarize(C_p_c_i * Dr(T_c(t, r)))...]
Ni = [scalarize(D_i_m .* Dr.(C_c_i))...]
Nidr = [scalarize(expand_Dr_D_im .* Dr.(C_c_i) + D_i_m .* Drr.(C_c_i))...]

# WGSParticle
eqs_Dim = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T_c(t, r), P))...]
eqs_ri = [scalarize(r_i .~ r_i_func(C_c_i, d_cat, θ, P, T_c(t, r)))...]
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) + r_i[i] for i in 1:5]

# WGSParticleExpanded
eqs_Hci = [scalarize(H_c_i .~ H_i_func(T_c(t, r)))...]
eqs_Cpci = [scalarize(C_p_c_i .~ C_p_i_func(T_c(t, r)))...]
eqs_Cpcat = [C_p_cat ~ C_p_cat_func(T_c(t, r))]
STEP_DE5_1 = sum([C_p_c_i[i] * C_c_i[i] for i in 1:5])
STEP_DE5_2 = sum([Dt(C_c_i[i]) * H_c_i[i] for i in 1:5])
STEP_DE5_3 = sum([H_c_i[i] * Ni[i] for i in 1:5])
STEP_DE5_4 = sum([H_c_i[i] * Nidr[i] for i in 1:5])
STEP_DE5_5 = sum([Ni[i] * Hidr[i] for i in 1:5])
STEP_DE5_6 = sum([H_c_i[i] * r_i[i] for i in 1:5])
DE5 = [(1 - θ) * ρ_cat * C_p_cat * Dt(T_c(t, r)) + θ * (STEP_DE5_1 * Dt(T_c(t, r)) + STEP_DE5_2) ~ (((λ_cat) / r) * Dr(T_c(t, r)) + λ_cat * Drr(T_c(t, r))) + θ * ( (2 / r) * STEP_DE5_3 + STEP_DE5_4 + STEP_DE5_5 - STEP_DE5_6)]
eqs = [eqs_Dim...; eqs_ri...; DE4...; eqs_Hci...; eqs_Cpci...; eqs_Cpcat; DE5]

# WGSParticle
ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2], C_c_3(0.0, r) ~ C_c_i_init[3], C_c_4(0.0, r) ~ C_c_i_init[4], C_c_5(0.0, r) ~ C_c_i_init[5]]
BCS2 = [Dr(C_c_1(t, 0.0)) ~ 0.0, Dr(C_c_2(t, 0.0)) ~ 0.0, Dr(C_c_3(t, 0.0)) ~ 0.0, Dr(C_c_4(t, 0.0)) ~ 0.0, Dr(C_c_5(t, 0.0)) ~ 0.0]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]

# WGSParticleExpanded
ICS_T_c = [T_c(0.0, r) ~ T_val]
BCS_Tc = [Dr(T_c(t, 0.0)) ~ 0.0]
STEP_BCS4_1 = sum([H_i[i] * k_c_i[i] * (C_c_i_rad[i] - C_i[i]) for i in 1:5])
STEP_BCS4_2 = sum([H_c_i[i] * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5])
BCS4 = [hf * (T_c(t, rad_cat) - T) + STEP_BCS4_1 ~ (-λ_cat) * Dr(T_c(t, rad_cat)) - STEP_BCS4_2]

bcs = [ICS_C_c_i...; BCS2...; BCS3...; ICS_T_c; BCS_Tc; BCS4]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain (time is in [h])
domains = [t ∈ Interval(0.0, 0.01),
    r ∈ Interval(0.0, rad_cat)]

# System
vars_ind = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), T_c(t, r)]
vars_dep = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5, C_p_c_1, C_p_c_2, C_p_c_3, C_p_c_4, C_p_c_5, H_c_1, H_c_2, H_c_3, H_c_4, H_c_5, C_p_cat]
vars = [vars_ind...; vars_dep...]
params_scal = [T => T_val, P => P_val, R => R_atmm3, θ => θ_val, τ => τ_val, d_cat => d_cat_val, ρ_cat => ρ_cat_val, λ_cat => λ_cat_val, hf => h_f_val]
params_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
params_vec_k_c_i = [k_c_i[i] => k_c_i_val[i] for i in 1:5]
params_vec_H_i = [H_i[i] => H_i_val[i] for i in 1:5]
params = [params_scal...; params_vec_C_i...; params_vec_k_c_i...; params_vec_H_i...;]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

WGS_pde.ivs

# Solving ODE
t0 = time()
sol = solve(prob, KenCarp47(), saveat = 0.0001, abstol = 1e-6, reltol = 1e-6)
t_1 = time() - t0
println("Time to solve ODE: ", t_1)
# # sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
# sols = sol[C_c_5(t, r)]

using DelimitedFiles

# Define T and P ranges
temp_range = [393.0; 483.0; 573.0;]
pres_range = [1.0; 2.0; 3.0;]
 
# Generate results
# make_results(500.1, 1.3, params, prob, 1e-8, 1e-12, 1e-12)

# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_results(temp_range[i], pres_range[j], params, prob)
#     end
# end

using Plots

# Generate plots
# make_plots(500.1, 1.3, [1; 11; 21], 1e-8, 1e-8)

# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_plots(temp_range[i], pres_range[j], [1; 11; 21], 0.001)
#     end
# end

# # Plotting 
# time = 0.01
# index_sol = Int(time/0.0001)
# solution = sols[1:index_sol, :]

# using Plots

# plot(solution)

# using DelimitedFiles
# folder = "WGS_particle_expanded/results_particle_expanded_lab"
# write_to_csv("C_c_1_lab.csv", sol[C_c_1(t, r)], folder)
# write_to_csv("C_c_2_lab.csv", sol[C_c_2(t, r)], folder)
# write_to_csv("C_c_3_lab.csv", sol[C_c_3(t, r)], folder)
# write_to_csv("C_c_4_lab.csv", sol[C_c_4(t, r)], folder)
# write_to_csv("C_c_5_lab.csv", sol[C_c_5(t, r)], folder)

