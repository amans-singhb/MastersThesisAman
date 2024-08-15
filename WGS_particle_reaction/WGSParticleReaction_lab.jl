##### Catalyst particle balance for WGS reactor model with reaction #####

using Pkg
Pkg.activate("WGS")
include("functions_WGSParticleReaction.jl")
include("functions_CTESN.jl")

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
T_val = 503.0 # [K] [param1]
P_val = 1.3 # [atm] [param6]

C_i_old = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]

y_0 = [0.208917; 0.0910445; 0.204129; 0.480558; 0.0153515;] # from paper
C_total = sum(C_i_old)
V = F_0/C_total

C_i_val_temp = (F_0 * y_0) / V # [mol/m^3]
C_i_val = [C_i_val_temp[1], C_i_val_temp[2], C_i_val_temp[3], 2.3 * C_i_val_temp[1], C_i_val_temp[5]] # [mol/m^3]

C_c_i_init = C_i_val
# zero_val = 1e-15 # isapprox(0, 1e-324) = true
# C_c_i_init = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]


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

# Domain postdiff (time is in [h])
domains = [t ∈ Interval(0.0, 2e-5),
    r ∈ Interval(0.0, rad_cat)]

# # Domain prediff (time is in [h])
# domains = [t ∈ Interval(0.0, 3e-3),
#     r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5]
prms_scal = [T => T_val, P => P_val, R => R_atmm3, θ => θ_val, τ => τ_val, d_cat => d_cat_val]
prms_vec_k_c_i = [k_c_i[i] => k_c_i_val[i] for i in 1:5]
prms_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
prms = [prms_scal...; prms_vec_k_c_i...; prms_vec_C_i...;]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, prms)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
t0_disc = time()
prob = discretize(WGS_pde, discretization)
t1_disc = time() - t0_disc
print("\nDiscretization time: ", t1_disc, "\n")

# Solving ODE
t0_sol = time()
sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
t1_sol = time() - t0_sol
print("\nSolution time: ", t1_sol, "\n")



# # Solving ODE
# t0_sol = time()

# # sol = solve(prob, FBDF(), abstol = 1e-6, reltol = 1e-6)
# # sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
# # sol = solve(prob, QNDF(), abstol = 1e-6, reltol = 1e-6)
# # sol = solve(prob, RadauIIA5(), abstol = 1e-6, reltol = 1e-6)
# # sol = solve(prob, Rodas5P(), abstol = 1e-6, reltol = 1e-6)

# t1_sol = time() - t0_sol
# print("\nSolution time: ", t1_sol, "\n")

# using DelimitedFiles

# timevec = []
# timevec = vec(readdlm("WGS_particle_reaction/solver_results/solvertimes_postdiff.csv", ',', '\n'))

# # time_string = "FBDF, $t1_sol"
# # time_string = "KenCarp47, $t1_sol"
# # time_string = "QNDF, $t1_sol"
# # time_string = "RadauIIA5, $t1_sol"
# # time_string = "Rodas5P, $t1_sol"

# push!(timevec, time_string)
# write_to_csv("solvertimes_postdiff.csv", timevec, "WGS_particle_reaction/solver_results")

# using JLD2
# using Plots

# sol_data = extractData(sol)

# # save("WGS_particle_reaction/solver_results/sol_FBDF_postdiff.jld2", "sol", sol_data)
# # test = load("WGS_particle_reaction/solver_results/sol_FBDF_postdiff.jld2", "sol")

# # save("WGS_particle_reaction/solver_results/sol_KenCarp47_postdiff.jld2", "sol", sol_data)
# # test = load("WGS_particle_reaction/solver_results/sol_KenCarp47_postdiff.jld2", "sol")

# # save("WGS_particle_reaction/solver_results/sol_QNDF_postdiff.jld2", "sol", sol_data)
# # test = load("WGS_particle_reaction/solver_results/sol_QNDF_postdiff.jld2", "sol")

# # save("WGS_particle_reaction/solver_results/sol_RadauIIA5_postdiff.jld2", "sol", sol_data)
# # test = load("WGS_particle_reaction/solver_results/sol_RadauIIA5_postdiff.jld2", "sol")

# # save("WGS_particle_reaction/solver_results/sol_Rodas5P_postdiff.jld2", "sol", sol_data)
# # test = load("WGS_particle_reaction/solver_results/sol_Rodas5P_postdiff.jld2", "sol")

# plot(sol.t*3600, test[21])
# plot!(sol.t*3600, sol_data[21])

#--------------------------------------------------------------------------------------------------------------------------#
# sols1 = sol[C_c_1(t, r)][:, 21]
# sols2 = sol[C_c_2(t, r)][:, 21]
# sols3 = sol[C_c_3(t, r)][:, 21]
# sols4 = sol[C_c_4(t, r)][:, 21]
# sols5 = sol[C_c_5(t, r)][:, 21]
# sol_t = sol.t * 3600

# plot(sol_t, sols1, label = "CO")
# plot!(sol_t, sols2, label = "CO2")
# plot!(sol_t, sols3, label = "H2")
# plot!(sol_t, sols4, label = "H2O")
# plot!(sol_t, sols5, label = "N2")

# fieldnames(typeof(sol))
# typeof(sol.original_sol)
# # sol = solve(prob, KenCarp47(), saveat = 1e-6, abstol = 1e-6, reltol = 1e-6)

# sol.t is decided by solver
# sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)

# using JLD2

# # sol.t length is 10
# sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
# lt = length(sol.t)
# sol.t
# soldata = extractData(sol)
# save("WGS_particle_reaction/sol_0/sol_0_t$lt.jld2", "sol", sol)

# # sol.t length is 21
# sol = solve(prob, KenCarp47(), saveat = 1e-6, abstol = 1e-6, reltol = 1e-6)
# lt = length(sol.t)
# soldata = extractData(sol)
# save("WGS_particle_reaction/sol_0/sol_0_t$lt.jld2", "sol", sol)

# # sol.t length is 101
# sol = solve(prob, KenCarp47(), saveat = 2e-7, abstol = 1e-6, reltol = 1e-6)
# lt = length(sol.t)
# soldata = extractData(sol)
# save("WGS_particle_reaction/sol_0/sol_0_t$lt.jld2", "sol", sol)

# # sol.t length is 1001
# sol = solve(prob, KenCarp47(), saveat = 2e-8, abstol = 1e-6, reltol = 1e-6)
# lt = length(sol.t)
# soldata = extractData(sol)
# save("WGS_particle_reaction/sol_0/sol_0_t$lt.jld2", "sol", sol)

# [soldata[i][1, 4] / soldata[i][1, 1] for i in 1:21]

# # sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)

# using DelimitedFiles

# Define T and P ranges
# temp_range = [393.0; 483.0; 573.0;]
# pres_range = [1.0; 2.0; 3.0;]
 
# Generate results
# make_results(500.1, 1.3, prms, prob)

# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_results(temp_range[i], pres_range[j], prms, prob)
#     end
# end

# using Plots

# Generate plots
# make_plots(500.1, 1.3, [1; 11; 21], 1e-8, 1e-8)

# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_plots(temp_range[i], pres_range[j], [1; 11; 21], 0.001)
#     end
# end

# ## Generate data for CTESN ##

# # Parameters
# y_p = [0.208917; 0.0910445; 0.204129; 0.480558; 0.0153515;] # from paper
# C_total = sum(C_i_old)
# V = F_0/C_total

# C_i_ml = (F_0 * y_p) / V # [mol/m^3]

# # C_1_ml = C_i_ml[1] # concentration of CO [mol/m^3]
# # T_ml = 503 # [K]
# # ratio_CO_H20 = 2.3

# # # Training parameters
# # p0 = [T_ml; ratio_CO_H20 * C_1_ml;]
# # p0_low = 0.9 * p0
# # p0_high = 1.1 * p0

# # Number of observations
# n_obs = 100

# # Training matrix

# # first time training (IMPORTANT: use line below parent_folder to save training parameters for future use)
# # p_train = (p0_low .+ (p0_high .- p0_low) .* rand(length(p0), n_obs))'

# # for regenerating data for same p_train (if needed)
# p_train = readdlm("WGS_particle_reaction/p_train.csv", ',', Float64, '\n')

# prms_ml = prms[2:end-5]
# prms_vec_C_i_ml = [C_i[i] => C_i_ml[i] for i in 1:3]
# prms_ml = [prms_ml...; prms_vec_C_i_ml...]

# using JLD2

# # folder_path_jld2 = "WGS_particle_reaction/ml_data_jld2_redone"
# parent_folder = "WGS_particle_reaction"
# # write_to_csv("p_train.csv", p_train, parent_folder)

# folder_time = "WGS_particle_reaction/training_data_time"
# times_vec = []
# t0_full = time()
# sim_full = []

# # Generating training data
# for i in 1:n_obs
#     t0 = time()

#     newprms_ml = prms_ml
#     newprms_ml = [T => p_train[i, 1]; prms_ml...; C_i[4] => p_train[i, 2]; C_i[5] => C_i_ml[5];]
#     newprob = remake(prob, p = newprms_ml)
#     newsol = solve(newprob, KenCarp47(), saveat = sol.t,abstol = 1e-6, reltol = 1e-6)
    
#     newsoldata = extractData(newsol)
#     push!(sim_full, newsoldata)
    
#     t1 = time() - t0
#     push!(times_vec, t1)
# end

# t1_full = time() - t0_full
# push!(times_vec, t1_full)
# print("\nTotal time: ", t1_full, "\n")
# write_to_csv("times_vec_t$lt.csv", times_vec, folder_time)

# size(sim_full)
# size(sim_full[1])
# size(sim_full[1][1])

# r_length = 21
# sim = []
# for i in 1:r_length
#     sim_r = []
#     for j in 1:size(p_train)[1]
#         push!(sim_r, sim_full[j][i])
#     end
#     push!(sim, sim_r)
# end

# size(sim)
# size(sim[1])
# size(sim[1][1])
# sim[1][1]
# test_ratio = [[sim[i][j][1, 4] / sim[i][j][1, 1] for j in 1:size(p_train)[1]] for i in 1:r_length]

# save("WGS_particle_reaction/sim/sim_t$lt.jld2", "sim", sim)
# test_sim = load("WGS_particle_reaction/sim/sim_t$lt.jld2", "sim")
# test_sim == sim
