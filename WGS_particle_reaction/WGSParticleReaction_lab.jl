##### Catalyst particle balance for WGS reactor model with reaction #####

using Pkg
Pkg.activate("WGS")
include("functions_WGSParticleReaction.jl")

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
# 1. T_c = T, and id kept constant, as we are only interested in the evolution of C_c_i
# 2. P_c = P, and is kept constant

### PDE system ###

## Parameters

# Dimensions
D_rct_val = 12.7e-3 # [m] 
rad_cat = 0.125e-3 # [m]
D_cat_val = rad_cat * 2 # [m]

# Catalyst properties
d_cat_val = 5904 # [kg/m^3] [param5]

# Coefficients
θ_val = 0.55 #[-] [param3]
τ_val = 5 #[-] [param4]
R_atmm3 = 8.2057e-2 # [m3 atm/kmol K] [param2]

# Inlet values
F_0 = 1e-5 # [mol/h] 
T_val = 500.1 # [K] [param1]
P_val = 1.3 # [atm] [param6]

C_i_val = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]
C_c_i_init = C_i_val * 0.75

# Mass transfer coefficients (bulk phase)
D_i_m_bulk = D_i_m_func(C_i_val, θ_val, τ_val, T_val, P_val) # [m^2/h]
k_c_i_val = k_c_i_func(T_val, P_val, R_atmm3, C_i_val, D_i_m_bulk, D_cat_val, D_rct_val, F_0)

using ModelingToolkit

## Parameters ##
@parameters t r T R θ τ d_cat P
@parameters k_c_i[1:5] C_i[1:5]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..)
@variables D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r) r_1(t, r) r_2(t, r) r_3(t, r) r_4(t, r) r_5(t, r)

C_c_i = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r)]
C_c_i_rad = [C_c_1(t, rad_cat), C_c_2(t, rad_cat), C_c_3(t, rad_cat), C_c_4(t, rad_cat), C_c_5(t, rad_cat)]
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]
r_i = [r_1, r_2, r_3, r_4, r_5]

## Equations and Differential Equations ##

# Differential of function D_i_m
Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T, P))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

using ModelingToolkit: scalarize

eqs_Dim = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T, P))...]
eqs_ri = [scalarize(r_i .~ r_i_func(C_c_i, d_cat, θ, P, T))...]
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) + r_i[i] for i in 1:5]


eqs = [eqs_Dim...; eqs_ri...; DE4...]

ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2], C_c_3(0.0, r) ~ C_c_i_init[3], C_c_4(0.0, r) ~ C_c_i_init[4], C_c_5(0.0, r) ~ C_c_i_init[5]]
BCS2 = [Dr(C_c_1(t, 0.0)) ~ 0.0, Dr(C_c_2(t, 0.0)) ~ 0.0, Dr(C_c_3(t, 0.0)) ~ 0.0, Dr(C_c_4(t, 0.0)) ~ 0.0, Dr(C_c_5(t, 0.0)) ~ 0.0]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]

bcs = [ICS_C_c_i...; BCS2...; BCS3...]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain (time is in [h])
domains = [t ∈ Interval(0.0, 0.00015),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5]
params_scal = [T => T_val, P => P_val, R => R_atmm3, θ => θ_val, τ => τ_val, d_cat => d_cat_val]
params_vec_k_c_i = [k_c_i[i] => k_c_i_val[i] for i in 1:5]
params_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
params = [params_scal...; params_vec_k_c_i...; params_vec_C_i...;]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

# Solving ODE
# sol = solve(prob, KenCarp47(), saveat = 0.000001, abstol = 1e-6, reltol = 1e-6)
# # sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
# sols = sol[C_c_1(t, r)]

using DelimitedFiles

# Define T and P ranges
# temp_range = [393.0; 483.0; 573.0;]
# pres_range = [1.0; 2.0; 3.0;]
 
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
# folder = "WGS_particle/results_particle_reaction_lab"
# write_to_csv("C_c_1_lab.csv", sol[C_c_1(t, r)], folder)
# write_to_csv("C_c_2_lab.csv", sol[C_c_2(t, r)], folder)
# write_to_csv("C_c_3_lab.csv", sol[C_c_3(t, r)], folder)
# write_to_csv("C_c_4_lab.csv", sol[C_c_4(t, r)], folder)
# write_to_csv("C_c_5_lab.csv", sol[C_c_5(t, r)], folder)


## Generate data for CTESN ##

# Parameters
y_p = [0.208917; 0.0910445; 0.204129; 0.480558; 0.0153515;] # from paper
C_total = sum(C_i_val)
V = F_0/C_total

C_i_ml = (F_0 * y_p) / V # [mol/m^3]

C_1_ml = C_i_ml[1] # concentration of CO [mol/m^3]
T_ml = 503 # [K]
ratio_CO_H20 = 2.3

p0 = [T_ml; ratio_CO_H20 * C_1_ml;]
p0_low = 0.9 * p0
p0_high = 1.1 * p0

# Number of observations
n_obs = 2

# training matrix
p_train = (p0_low .+ (p0_high .- p0_low) .* rand(length(p0), n_obs))'

params_ml = params[2:end-2]
params_vec_C_i = [C_i[i] => C_i_ml[i] for i in 1:5]

parent_folder = "WGS_particle_reaction/ml_data"

for i in 1:n_obs
    newparams_ml = params_ml
    newparams_ml = [T => p_train[i, 1]; params_ml...; C_i[4] => p_train[i, 2]; C_i[5] => C_i_ml[5];]
    newprob = remake(prob, p = newparams_ml)
    newsol = solve(newprob, KenCarp47(), saveat = 0.000001, abstol = 1e-6, reltol = 1e-6)

    string_param = string(p_train[i, 1]) * "K_" * string(p_train[i, 1]) * "molm3"
    folder = parent_folder * "/" * string_param

    string_cc1 = "C_c_1_" * string_param * "_lab.csv"
    string_cc2 = "C_c_2_" * string_param * "_lab.csv"
    string_cc3 = "C_c_3_" * string_param * "_lab.csv"
    string_cc4 = "C_c_4_" * string_param * "_lab.csv"
    string_cc5 = "C_c_5_" * string_param * "_lab.csv"

    write_to_csv(string_cc1, newsol[C_c_1(t, r)], folder)
    write_to_csv(string_cc2, newsol[C_c_2(t, r)], folder)
    write_to_csv(string_cc3, newsol[C_c_3(t, r)], folder)
    write_to_csv(string_cc4, newsol[C_c_4(t, r)], folder)
    write_to_csv(string_cc5, newsol[C_c_5(t, r)], folder)
end