##### Catalyst particle balance for WGS reactor model without reaction #####

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

C_i_val = (F_0 * y_0) / V

zero_val = 1e-15 # isapprox(0, 1e-324) = true
C_c_i_init = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]

# Mass transfer coefficients (bulk phase)
D_i_m_bulk = D_i_m_func(C_i_val, θ_val, τ_val, T_val, P_val) # [m^2/h]
k_c_i_val = k_c_i_func(T_val, P_val, R_atmm3, C_i_val, D_i_m_bulk, D_cat_val, D_rct_val, F_0)

using ModelingToolkit

## Parameters ##
@parameters t r T P R θ τ d_cat
@parameters C_i[1:5] k_c_i[1:5]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_1(..) C_c_2(..) C_c_3(..) C_c_4(..) C_c_5(..)
# dCc_1(..) dCc_2(..) dCc_3(..) dCc_4(..) dCC_5(..)
@variables D_1_m(t, r) D_2_m(t, r) D_3_m(t, r) D_4_m(t, r) D_5_m(t, r)

# dCc_i = [dCc_1(t, r), dCc_2(t, r), dCc_3(t, r), dCc_4(t, r), dCC_5(t, r)]
C_c_i = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r)]
C_c_i_rad = [C_c_1(t, rad_cat), C_c_2(t, rad_cat), C_c_3(t, rad_cat), C_c_4(t, rad_cat), C_c_5(t, rad_cat)]
D_i_m = [D_1_m, D_2_m, D_3_m, D_4_m, D_5_m]

## Equations and Differential Equations ##

# Differential of function D_i_m
Dr_D_im = Dr.(D_i_m_func(C_c_i, θ, τ, T, P))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

using ModelingToolkit: scalarize

# eqs_deriv = [dCc_i[i] ~ Dt(C_c_i[i]) for i in 1:5]
eqs_Dim = [scalarize(D_i_m .~ D_i_m_func(C_c_i, θ, τ, T, P))...]
DE4 = [Dt(C_c_i[i]) ~ (((2 * D_i_m[i]) / r) + expand_Dr_D_im[i]) * Dr(C_c_i[i]) + D_i_m[i] * Drr(C_c_i[i]) for i in 1:5]

eqs = [eqs_Dim...; DE4...]

ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2], C_c_3(0.0, r) ~ C_c_i_init[3], C_c_4(0.0, r) ~ C_c_i_init[4], C_c_5(0.0, r) ~ C_c_i_init[5]]
BCS2 = [Dr(C_c_1(t, 0.0)) ~ 0.0, Dr(C_c_2(t, 0.0)) ~ 0.0, Dr(C_c_3(t, 0.0)) ~ 0.0, Dr(C_c_4(t, 0.0)) ~ 0.0, Dr(C_c_5(t, 0.0)) ~ 0.0]
BCS3 = [k_c_i[i] * (C_c_i_rad[i] - C_i[i]) ~ (-1) * D_i_m[i] * Dr(C_c_i_rad[i]) for i in 1:5]

bcs = [ICS_C_c_i...; BCS2...; BCS3...]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain (time is in [h])
domains = [t ∈ Interval(0.0, 1e-2),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m] # dCc_1(t, r), dCc_2(t, r), dCc_3(t, r), dCc_4(t, r), dCC_5(t, r)]
prms_scal = [T => T_val, P => P_val, R => R_atmm3, θ => θ_val, τ => τ_val, d_cat => d_cat_val]
prms_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
prms_vec_k_c_i = [k_c_i[i] => k_c_i_val[i] for i in 1:5]
prms = [prms_scal...; prms_vec_C_i...; prms_vec_k_c_i...;]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, prms)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
# t0_disc = time()
prob = discretize(WGS_pde, discretization)
# t1_disc = time() - t0_disc
# println("Discretization time: ", t1_disc)

using BenchmarkTools

t0_solve = time()
sol = solve(prob, Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
t1_solve = time() - t0_solve
println("Solving time: ", t1_solve)

# # @btime
# @btime sol = solve(prob, FBDF(), abstol = 1e-6, reltol = 1e-6)
# @btime sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
# @btime sol = solve(prob, QNDF(), abstol = 1e-6, reltol = 1e-6)
# @btime sol = solve(prob, RadauIIA5(), abstol = 1e-6, reltol = 1e-6)
# @btime sol = solve(prob, Rodas5P(), abstol = 1e-6, reltol = 1e-6)

# @btime sol = solve(prob, Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)

# # @benchmark
# @benchmark sol = solve(prob, FBDF(), abstol = 1e-6, reltol = 1e-6)
# @benchmark sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
# @benchmark sol = solve(prob, QNDF(), abstol = 1e-6, reltol = 1e-6)
# @benchmark sol = solve(prob, RadauIIA5(), abstol = 1e-6, reltol = 1e-6)
# @benchmark sol = solve(prob, Rodas5P(), abstol = 1e-6, reltol = 1e-6)

# @benchmark sol = solve(prob, Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)


#--------------------------------------------------------------------------------------------------------------#

# sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
# sols = sol[C_c_1(t, r)]
# sol_t = sol.t * 3600

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

# using DelimitedFiles

# Define T and P rangues
# temp_range = [393.0; 483.0; 573.0;]
# pres_range = [1.0; 2.0; 3.0;]
 
# Generate results
# make_results(temp_range[3], pres_range[3], prms, prob)

# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_results(temp_range[i], pres_range[j], prms, prob)
#     end
# end

# using Plots

# Generate plots
# for i in eachindex(temp_range)
#     for j in eachindex(pres_range)
#         make_plots(temp_range[i], pres_range[j], [1; 11; 21], 0.001)
#     end
# end

# # Plotting 
# time = 0.01
# index_sol = Int(time/0.0001)
# solution = sols[:, 1]

# using Plots

# plot(sol_t, solution)

# folder = "WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm"
# write_to_csv("C_c_1_573.0K_3.0atm_lab.csv", sol[C_c_1(t, r)], folder)
# write_to_csv("C_c_2_573.0K_3.0atm_lab.csv", sol[C_c_2(t, r)], folder)
# write_to_csv("C_c_3_573.0K_3.0atm_lab.csv", sol[C_c_3(t, r)], folder)
# write_to_csv("C_c_4_573.0K_3.0atm_lab.csv", sol[C_c_4(t, r)], folder)
# write_to_csv("C_c_5_573.0K_3.0atm_lab.csv", sol[C_c_5(t, r)], folder)

# C_c_11 = readdlm("WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm/C_c_1_573.0K_3.0atm_lab.csv", ',', Float64, '\n')
# C_c_21 = readdlm("WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm/C_c_2_573.0K_3.0atm_lab.csv", ',', Float64, '\n')
# C_c_31 = readdlm("WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm/C_c_3_573.0K_3.0atm_lab.csv", ',', Float64, '\n')
# C_c_41 = readdlm("WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm/C_c_4_573.0K_3.0atm_lab.csv", ',', Float64, '\n')
# C_c_51 = readdlm("WGS_particle/results_particle_lab_parameters/param573.0K_3.0atm/C_c_5_573.0K_3.0atm_lab.csv", ',', Float64, '\n')
