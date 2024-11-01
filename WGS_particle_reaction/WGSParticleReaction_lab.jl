##### Catalyst particle balance for WGS reactor model with reaction #####

using Pkg
Pkg.activate("WGS")
Pkg.status()
include("functions_WGSParticleReaction.jl")
include("ALAMO_Surrogates.jl")

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

## Constant parameters

# Dimensions reactor
D_rct = 12.7e-3 # [m] 
L_rct = 4.8e-3 # [m]
V_rct = π * (D_rct/2)^2 * L_rct # [m^3]
L_rct2 = 304.8e-3 # [m]
V_rct2 = π * (D_rct/2)^2 * L_rct2 # [m^3]

# Dimensions particle
rad_cat = 0.125e-3 # [m]
D_cat = rad_cat * 2 # [m]

# Catalyst properties
d_cat = 5904 # [kg/m^3]

# Coefficients
θ = 0.55 #[-]
τ = 5 #[-]
R = 8.2057e-2 # [m3 atm/kmol K]

# Industrial reactor parameter from paper ------------------------------------- #
F_ind = 9199 * 1000 # [mol/h]
D_ind = 4.4 # [m]
L_ind = 14 # [m]

V_ind = π * (D_ind/2)^2 * L_ind # [m^3]
P_ind = 54.28 # [atm]
T_ind = 503.0 # [K]

q_ind = (F_ind * R * 1e-3 * T_ind) / P_ind # [m^3/h]
restime_ind = V_ind / q_ind # [h]
restime_ind_s = restime_ind * 3600 # [s]
# ----------------------------------------------------------------------------- #

## Parameters

# Inlet values
T_val = 503.0 # [K] [param1]
P_val = 1.3 # [atm] [param2]

# Temp vals F ----------------------------------------------------------------- #
y_0 = [0.208917; 0.0910445; 0.204129; 0.480558; 0.0153515;] # from paper
volume_ratio = V_rct / V_ind
vol_ratio2 = V_rct2 / V_ind
f_temp2 = F_ind * vol_ratio2
f_temp = F_ind * volume_ratio
# ----------------------------------------------------------------------------- #
F_val = f_temp * (P_val/ P_ind) * (T_ind / T_val) # [mol/h] [param3]

# Temp vals C ----------------------------------------------------------------- #
C_total = P_val / (R * 1e-3 * T_val) # [mol/m^3]
C_i_val_temp = C_total * y_0 # [mol/m^3]
C_i_val_temp[4]/C_i_val_temp[1] # 2.3
# ----------------------------------------------------------------------------- #
C_i_val = [C_i_val_temp[1], C_i_val_temp[2], C_i_val_temp[3], 2.3 * C_i_val_temp[1], C_i_val_temp[5]] # [mol/m^3]

# # postdiff:
# C_c_i_init_val = C_i_val

# prediff:
zero_val = 1e-15 # isapprox(0, 1e-324) = true
C_c_i_init_val = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]

using ModelingToolkit

## Parameters ##
@parameters t r T P F
@parameters C_i[1:5] C_c_i_init[1:5]

# Mass transfer coefficients (bulk phase)
D_i_m_bulk = D_i_m_func(C_i, θ, τ, T, P) # [m^2/h]
k_c_i = k_c_i_func(T, P, R, C_i, D_i_m_bulk, D_cat, D_rct, F)

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

# # Domain postdiff (time is in [h])
# domains = [t ∈ Interval(0.0, 2e-5),
#     r ∈ Interval(0.0, rad_cat)]

# Domain prediff (time is in [h])
domains = [t ∈ Interval(0.0, 5e-4),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r), C_c_3(t, r), C_c_4(t, r), C_c_5(t, r), D_1_m, D_2_m, D_3_m, D_4_m, D_5_m, r_1, r_2, r_3, r_4, r_5]
prms_scal = [T => T_val, P => P_val, F => F_val]
prms_vec_C_i = [C_i[i] => C_i_val[i] for i in 1:5]
prms_vec_C_c_i_init = [C_c_i_init[i] => C_c_i_init_val[i] for i in 1:5]
prms = [prms_scal...; prms_vec_C_i...; prms_vec_C_c_i_init...]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, prms)

# Discretization
dr = rad_cat/20
order = 2
discretization = MOLFiniteDifference([r => dr], t, order=order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)

# Solving ODE
sol = solve(prob, FBDF(), saveat = 1e-6, abstol = 1e-6, reltol = 1e-6)

using Plots

sols1 = sol[C_c_1(t, r)][:, 1]
sols2 = sol[C_c_2(t, r)][:, 1]
sols3 = sol[C_c_3(t, r)][:, 1]
sols4 = sol[C_c_4(t, r)][:, 1]
sols5 = sol[C_c_5(t, r)][:, 1]
solt = sol.t * 3600
sol_t = solt

plot(sol_t, sols1, label = "C_c_1", dpi = 1000)
plot!(sol_t, sols2, label = "C_c_2")
plot!(sol_t, sols3, label = "C_c_3")
plot!(sol_t, sols4, label = "C_c_4")
plot!(sol_t, sols5, label = "C_c_5")
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
title!("True concentration profiles of C_c_i for case from paper", titlefontsize=10, titlefontcolor=:black)
savefig("WGS_particle_reaction/figures_rslt/result_true_case1.png")

# ------------------------------------------------------------------------------------------------------------------------------------------- #
prms_int = readdlm("WGS_particle_reaction/csv_files/samples_LHS_init.csv", ',', Float64, '\n')
prms_beyond = readdlm("WGS_particle_reaction/csv_files/samples_LHS_beyond.csv", ',', Float64, '\n')
prms_random = readdlm("WGS_particle_reaction/csv_files/samples_LHS_random.csv", ',', Float64, '\n')

cci_bnd = readdlm("WGS_particle_reaction/csv_files/cci_t0_beyond.csv", ',', Float64, '\n')

cci_beyond_1 = repeat(cci_bnd[:, 1], inner = 100)
cci_beyond_2 = repeat(cci_bnd[:, 2], inner = 100)
cci_beyond_3 = repeat(cci_bnd[:, 3], inner = 100)
cci_beyond_4 = repeat(cci_bnd[:, 4], inner = 100)
cci_beyond_5 = repeat(cci_bnd[:, 5], inner = 100)
cci_beyond = [cci_beyond_1 cci_beyond_2 cci_beyond_3 cci_beyond_4 cci_beyond_5]

cci_int_1 = repeat([C_c_i_init_val[1]], 100)
cci_int_2 = repeat([C_c_i_init_val[2]], 100)
cci_int_3 = repeat([C_c_i_init_val[3]], 100)
cci_int_4 = repeat([C_c_i_init_val[4]], 100)
cci_int_5 = repeat([C_c_i_init_val[5]], 100)
cci_init = [cci_int_1 cci_int_2 cci_int_3 cci_int_4 cci_int_5]

input_init = [prms_int cci_init]
input_beyond = [prms_beyond cci_beyond]
input = [input_init; input_beyond; prms_random]
# ------------------------------------------------------------------------------------------------------------------------------------------- #

# Using surrogate model to predict the evolution of C_c_i with random parameters

using Random

T_val_r = rand(393.0:573.0)

C_1_r = rand(Float64)
C_2_r = rand(Float64)
C_3_r = rand(Float64)
C_4_r = rand(Float64)
C_5_r = rand(Float64)
C_i_r = [C_1_r, C_2_r, C_3_r, C_4_r, C_5_r]
C_i_r = (C_i_r ./ sum(C_i_r)) * C_total

C_c_1_init_r = rand(Float64)
C_c_2_init_r = rand(Float64)
C_c_3_init_r = rand(Float64)
C_c_4_init_r = rand(Float64)
C_c_5_init_r = rand(Float64)
C_c_i_init_r = [C_c_1_init_r, C_c_2_init_r, C_c_3_init_r, C_c_4_init_r, C_c_5_init_r]
C_c_i_init_r = (C_c_i_init_r ./ sum(C_c_i_init_r)) * C_total

for i in 1:11100
    if input[i, 1] == T_val_r && input[i, 2] == C_1_r && input[i, 3] == C_2_r && input[i, 4] == C_3_r && input[i, 5] == C_4_r && input[i, 6] == C_5_r && input[i, 7] == C_c_1_init_r && input[i, 8] == C_c_2_init_r && input[i, 9] == C_c_3_init_r && input[i, 10] == C_c_4_init_r && input[i, 11] == C_c_5_init_r
        println(i)
    end
end

save_vec = [T_val_r, C_1_r, C_2_r, C_3_r, C_4_r, C_5_r, C_c_1_init_r, C_c_2_init_r, C_c_3_init_r, C_c_4_init_r, C_c_5_init_r]
writedlm("WGS_particle_reaction/csv_files/random_params.csv", save_vec, ',')

prms_scal_r = [T => T_val_r, P => P_val, F => F_val]
prms_C_i_r = [C_i[i] => C_i_r[i] for i in 1:5]
prms_C_c_i_r = [C_c_i_init[i] => C_c_i_init_r[i] for i in 1:5]
prms_r = [prms_scal_r...; prms_C_i_r...; prms_C_c_i_r...]
@named WGS_pde_r = PDESystem(eqs, bcs, domains, [t, r], vars, prms_r)

prob_r = discretize(WGS_pde_r, discretization)
sol_true_r = solve(prob_r, FBDF(), saveat = 1e-6 ,abstol = 1e-6, reltol = 1e-6)

sols1_r = sol_true_r[C_c_1(t, r)][:, 1]
sols2_r = sol_true_r[C_c_2(t, r)][:, 1]
sols3_r = sol_true_r[C_c_3(t, r)][:, 1]
sols4_r = sol_true_r[C_c_4(t, r)][:, 1]
sols5_r = sol_true_r[C_c_5(t, r)][:, 1]

plot(sol_t, sols1_r, label = "C_c_1", dpi = 1000)
plot!(sol_t, sols2_r, label = "C_c_2")
plot!(sol_t, sols3_r, label = "C_c_3")
plot!(sol_t, sols4_r, label = "C_c_4")
plot!(sol_t, sols5_r, label = "C_c_5")
title!("True concentration profiles of C_c_i for random case", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
savefig("WGS_particle_reaction/figures_rslt/result_true_case_rand.png")

# ------------------------------------------------------------------------------------------------------------------------------------------- #
# Using surrogate model to preedict the evolution of C_c_i

using DifferentialEquations

function WGS_surrogate!(du, u, p, t)
    du[1] = dcc1dt_surrogate(p[1], p[2], p[3], p[4], p[5], p[6], u[1], u[2], u[3], u[4], u[5])
    du[2] = dcc2dt_surrogate(p[1], p[2], p[3], p[4], p[5], p[6], u[1], u[2], u[3], u[4], u[5])
    du[3] = dcc3dt_surrogate(p[1], p[2], p[3], p[4], p[5], p[6], u[1], u[2], u[3], u[4], u[5])
    du[4] = dcc4dt_surrogate(p[1], p[2], p[3], p[4], p[5], p[6], u[1], u[2], u[3], u[4], u[5])
    du[5] = dcc5dt_surrogate(p[1], p[2], p[3], p[4], p[5], p[6], u[1], u[2], u[3], u[4], u[5])
end

tspan = (0.0, 5e-4)

## Case from paper (not in the training data) ##

T_val = 503.0

C_1 = 6.580115441743921
C_2 = 2.8675661642463486
C_3 = 6.429311090087188
C_4 = 15.134265516011016
C_5 = 0.4835156650915522

C_1 == C_i_val[1] # true
C_2 == C_i_val[2] # true
C_3 == C_i_val[3] # true
C_4 == C_i_val[4] # true
C_5 == C_i_val[5] # true

for i in 1:11100
    if input[i, 1] == T_val && input[i, 2] == C_1 && input[i, 3] == C_2 && input[i, 4] == C_3 && input[i, 5] == C_4 && input[i, 6] == C_5 && input[i, 7] == C_c_i_init_val[1] && input[i, 8] == C_c_i_init_val[2] && input[i, 9] == C_c_i_init_val[3] && input[i, 10] == C_c_i_init_val[4] && input[i, 11] == C_c_i_init_val[5]
        println(i)
    end
end

u0 = [1.0e-15, 1.0e-15, 1.0e-15, 1.0e-15, 31.49631404693692]
p = [T_val, C_1, C_2, C_3, C_4, C_5]
prob = ODEProblem(WGS_surrogate!, u0, tspan, p)
sol_pred = solve(prob, Tsit5(), saveat = 1e-6)

sol_pred_m = (sol_pred[1:end, 1:end])'

plot(sol_t, sol_pred_m, label=["C_c_1" "C_c_2" "C_c_3" "C_c_4" "C_c_5"], dpi = 1000)
title!("Predicted concentration profiles of C_c_i for case from paper", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
savefig("WGS_particle_reaction/figures_rslt/result_pred_case1.png")

error1 = abs.(sol_pred[1, :] - sols1)
error2 = abs.(sol_pred[2, :] - sols2)
error3 = abs.(sol_pred[3, :] - sols3)
error4 = abs.(sol_pred[4, :] - sols4)
error5 = abs.(sol_pred[5, :] - sols5)

plot(sol_t, error1, label="Error C_c_1", dpi = 1000)
plot!(sol_t, error2, label="Error C_c_2")
plot!(sol_t, error3, label="Error C_c_3")
plot!(sol_t, error4, label="Error C_c_4")
plot!(sol_t, error5, label="Error C_c_5")
title!("Error in prediction of C_c_i for case from paper", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Error [mol/m^3]")
savefig("WGS_particle_reaction/figures_rslt/error_case1.png")

## Random case ##

u0_r = C_c_i_init_r
p_r = [T_val_r, C_i_r[1], C_i_r[2], C_i_r[3], C_i_r[4], C_i_r[5]]
prob_r = ODEProblem(WGS_surrogate!, u0_r, tspan, p_r)
sol_pred_r = solve(prob_r, Tsit5(), saveat = 1e-6)

sol_pred_r_m = (sol_pred_r[1:end, 1:end])'

plot(sol_t, sol_pred_r_m, label=["C_c_1" "C_c_2" "C_c_3" "C_c_4" "C_c_5"], dpi = 1000)
title!("Predicted concentration profiles of C_c_i for random case", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
savefig("WGS_particle_reaction/figures_rslt/result_pred_case_rand.png")

error1_r = abs.(sol_pred_r[1, :] - sols1_r)
error2_r = abs.(sol_pred_r[2, :] - sols2_r)
error3_r = abs.(sol_pred_r[3, :] - sols3_r)
error4_r = abs.(sol_pred_r[4, :] - sols4_r)
error5_r = abs.(sol_pred_r[5, :] - sols5_r)

plot(sol_t, error1_r, label="Error C_c_1", dpi = 1000)
plot!(sol_t, error2_r, label="Error C_c_2")
plot!(sol_t, error3_r, label="Error C_c_3")
plot!(sol_t, error4_r, label="Error C_c_4")
plot!(sol_t, error5_r, label="Error C_c_5")
title!("Error in prediction of C_c_i for random case", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Error [mol/m^3]")
savefig("WGS_particle_reaction/figures_rslt/error_case_rand.png")

# ------------------------------------------------------------------------------------------------------------------------------------------- #
# Using surrogate model to predict points
using DelimitedFiles

# Load test data and training data, and divide it into input and output
test_data = readdlm("WGS_particle_reaction/csv_files/test_data_m2.csv", ',', Float64, '\n')
train_data = readdlm("WGS_particle_reaction/csv_files/all_data_m2.csv", ',', Float64, '\n')

test_input = test_data[:, 1:11]
true_output = test_data[:, 12:end]

train_input = train_data[:, 1:11]
train_output = train_data[:, 12:end]

# Predicting the output using the surrogate models
pred_output_dcc1 = dcc1dt_surrogate3.(test_input[:, 1], test_input[:, 2], test_input[:, 3], test_input[:, 4], test_input[:, 5], test_input[:, 6], test_input[:, 7], test_input[:, 8], test_input[:, 9], test_input[:, 10], test_input[:, 11]) 
pred_output_dcc2 = dcc2dt_surrogate3.(test_input[:, 1], test_input[:, 2], test_input[:, 3], test_input[:, 4], test_input[:, 5], test_input[:, 6], test_input[:, 7], test_input[:, 8], test_input[:, 9], test_input[:, 10], test_input[:, 11])
pred_output_dcc3 = dcc3dt_surrogate3.(test_input[:, 1], test_input[:, 2], test_input[:, 3], test_input[:, 4], test_input[:, 5], test_input[:, 6], test_input[:, 7], test_input[:, 8], test_input[:, 9], test_input[:, 10], test_input[:, 11])
pred_output_dcc4 = dcc4dt_surrogate3.(test_input[:, 1], test_input[:, 2], test_input[:, 3], test_input[:, 4], test_input[:, 5], test_input[:, 6], test_input[:, 7], test_input[:, 8], test_input[:, 9], test_input[:, 10], test_input[:, 11])
pred_output_dcc5 = dcc5dt_surrogate3.(test_input[:, 1], test_input[:, 2], test_input[:, 3], test_input[:, 4], test_input[:, 5], test_input[:, 6], test_input[:, 7], test_input[:, 8], test_input[:, 9], test_input[:, 10], test_input[:, 11])
pred_output = [pred_output_dcc1 pred_output_dcc2 pred_output_dcc3 pred_output_dcc4 pred_output_dcc5]

train_output_dcc1 = dcc1dt_surrogate3.(train_input[:, 1], train_input[:, 2], train_input[:, 3], train_input[:, 4], train_input[:, 5], train_input[:, 6], train_input[:, 7], train_input[:, 8], train_input[:, 9], train_input[:, 10], train_input[:, 11])
train_output_dcc2 = dcc2dt_surrogate3.(train_input[:, 1], train_input[:, 2], train_input[:, 3], train_input[:, 4], train_input[:, 5], train_input[:, 6], train_input[:, 7], train_input[:, 8], train_input[:, 9], train_input[:, 10], train_input[:, 11])
train_output_dcc3 = dcc3dt_surrogate3.(train_input[:, 1], train_input[:, 2], train_input[:, 3], train_input[:, 4], train_input[:, 5], train_input[:, 6], train_input[:, 7], train_input[:, 8], train_input[:, 9], train_input[:, 10], train_input[:, 11])
train_output_dcc4 = dcc4dt_surrogate3.(train_input[:, 1], train_input[:, 2], train_input[:, 3], train_input[:, 4], train_input[:, 5], train_input[:, 6], train_input[:, 7], train_input[:, 8], train_input[:, 9], train_input[:, 10], train_input[:, 11])
train_output_dcc5 = dcc5dt_surrogate3.(train_input[:, 1], train_input[:, 2], train_input[:, 3], train_input[:, 4], train_input[:, 5], train_input[:, 6], train_input[:, 7], train_input[:, 8], train_input[:, 9], train_input[:, 10], train_input[:, 11])
pred_train_output = [train_output_dcc1 train_output_dcc2 train_output_dcc3 train_output_dcc4 train_output_dcc5]

# Quality metrics dCC1/dt
mse_dcc1 = MSE(true_output[:, 1], pred_output_dcc1)
rmse_dcc1 = RMSE(true_output[:, 1], pred_output_dcc1)
r_sqr_dcc1 = R_squared(true_output[:, 1], pred_output_dcc1)

mse_dcc1_train = MSE(train_output[:, 1], train_output_dcc1)
rmse_dcc1_train = RMSE(train_output[:, 1], train_output_dcc1)
r_sqr_dcc1_train = R_squared(train_output[:, 1], train_output_dcc1)

## dCC1/dt
plot(true_output[:, 1], pred_output_dcc1, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc1, digits = 3))", dpi = 600)
plot!(true_output[:, 1], true_output[:, 1], ls = :dash, label = false)
title!("Parity plot for dC_c_1/dt using test data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc1.png")

plot(train_output[:, 1], train_output_dcc1, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc1_train, digits = 3))", dpi = 600)
plot!(train_output[:, 1], train_output[:, 1], ls = :dash, label = false)
title!("Parity plot for dC_c_1/dt using training data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc1_train.png")

# Quality metrics dCC2/dt
mse_dcc2 = MSE(true_output[:, 2], pred_output_dcc2)
rmse_dcc2 = RMSE(true_output[:, 2], pred_output_dcc2)
r_sqr_dcc2 = R_squared(true_output[:, 2], pred_output_dcc2)

mse_dcc2_train = MSE(train_output[:, 2], train_output_dcc2)
rmse_dcc2_train = RMSE(train_output[:, 2], train_output_dcc2)
r_sqr_dcc2_train = R_squared(train_output[:, 2], train_output_dcc2)

# dCC2/dt
plot(true_output[:, 2], pred_output_dcc2, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc2, digits = 3))", dpi = 600)
plot!(true_output[:, 2], true_output[:, 2], ls = :dash, label = false)
title!("Parity plot for dC_c_2/dt using test data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc2.png")

plot(train_output[:, 2], train_output_dcc2, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc2_train, digits = 3))", dpi = 600)
plot!(train_output[:, 2], train_output[:, 2], ls = :dash, label = false)
title!("Parity plot for dC_c_2/dt using training data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc2_train.png")

# Quality metrics dCC3/dt
mse_dcc3 = MSE(true_output[:, 3], pred_output_dcc3)
rmse_dcc3 = RMSE(true_output[:, 3], pred_output_dcc3)
r_sqr_dcc3 = R_squared(true_output[:, 3], pred_output_dcc3)

mse_dcc3_train = MSE(train_output[:, 3], train_output_dcc3)
rmse_dcc3_train = RMSE(train_output[:, 3], train_output_dcc3)
r_sqr_dcc3_train = R_squared(train_output[:, 3], train_output_dcc3)

# dCC3/dt
plot(true_output[:, 3], pred_output_dcc3, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc3, digits = 3))", dpi = 600)
plot!(true_output[:, 3], true_output[:, 3], ls = :dash, label = false)
title!("Parity plot for dC_c_3/dt using test data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc3.png")

plot(train_output[:, 3], train_output_dcc3, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc3_train, digits = 3))", dpi = 600)
plot!(train_output[:, 3], train_output[:, 3], ls = :dash, label = false)
title!("Parity plot for dC_c_3/dt using training data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc3_train.png")

# Quality metrics dCC4/dt
mse_dcc4 = MSE(true_output[:, 4], pred_output_dcc4)
rmse_dcc4 = RMSE(true_output[:, 4], pred_output_dcc4)
r_sqr_dcc4 = R_squared(true_output[:, 4], pred_output_dcc4)

mse_dcc4_train = MSE(train_output[:, 4], train_output_dcc4)
rmse_dcc4_train = RMSE(train_output[:, 4], train_output_dcc4)
r_sqr_dcc4_train = R_squared(train_output[:, 4], train_output_dcc4)

# dCC4/dt
plot(true_output[:, 4], pred_output_dcc4, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc4, digits = 3))", dpi = 600)
plot!(true_output[:, 4], true_output[:, 4], ls = :dash, label = false)
title!("Parity plot for dC_c_4/dt using test data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc4.png")

plot(train_output[:, 4], train_output_dcc4, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc4_train, digits = 3))", dpi = 600)
plot!(train_output[:, 4], train_output[:, 4], ls = :dash, label = false)
title!("Parity plot for dC_c_4/dt using training data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc4_train.png")

# Quality metrics dCC5/dt
mse_dcc5 = MSE(true_output[:, 5], pred_output_dcc5)
rmse_dcc5 = RMSE(true_output[:, 5], pred_output_dcc5)
r_sqr_dcc5 = R_squared(true_output[:, 5], pred_output_dcc5)

mse_dcc5_train = MSE(train_output[:, 5], train_output_dcc5)
rmse_dcc5_train = RMSE(train_output[:, 5], train_output_dcc5)
r_sqr_dcc5_train = R_squared(train_output[:, 5], train_output_dcc5)

# dCC5/dt
plot(true_output[:, 5], pred_output_dcc5, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc5, digits = 3))", dpi = 600)
plot!(true_output[:, 5], true_output[:, 5], ls = :dash, label = false)
title!("Parity plot for dC_c_5/dt using test data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc5.png")

plot(train_output[:, 5], train_output_dcc5, seriestype = :scatter, label = "R² = $(round(r_sqr_dcc5_train, digits = 3))", dpi = 600)
plot!(train_output[:, 5], train_output[:, 5], ls = :dash, label = false)
title!("Parity plot for dC_c_5/dt using training data")
xlabel!("True values")
ylabel!("Predicted values")
savefig("WGS_particle_reaction/parity_plots_pred/parity_dcc5_train.png")

# r_sqr_train = round.([r_sqr_dcc1_train, r_sqr_dcc2_train, r_sqr_dcc3_train, r_sqr_dcc4_train, r_sqr_dcc5_train], digits = 3)
# r_sqr_test = round.([r_sqr_dcc1, r_sqr_dcc2, r_sqr_dcc3, r_sqr_dcc4, r_sqr_dcc5], digits = 3)

# mse_train = [mse_dcc1_train, mse_dcc2_train, mse_dcc3_train, mse_dcc4_train, mse_dcc5_train]
# mse_test = [mse_dcc1, mse_dcc2, mse_dcc3, mse_dcc4, mse_dcc5]

# rmse_train = [rmse_dcc1_train, rmse_dcc2_train, rmse_dcc3_train, rmse_dcc4_train, rmse_dcc5_train]
# rmse_test = [rmse_dcc1, rmse_dcc2, rmse_dcc3, rmse_dcc4, rmse_dcc5]

# q_train = [r_sqr_train mse_train rmse_train]
# q_test = [r_sqr_test mse_test rmse_test]

# write_to_csv("quality_metrics_train.csv", q_train, "WGS_particle_reaction/csv_files/")
# write_to_csv("quality_metrics_test.csv", q_test, "WGS_particle_reaction/csv_files/")

# # quality_metrics = transpose(hcat(r_sqr_train, r_sqr_test, mse_train, mse_test, rmse_train, rmse_test))

# # q = readdlm("WGS_particle_reaction/csv_files/quality_metrics_surrogate_models.csv", ',', Float64, '\n')

#-------------------------------------------------------------------------------------------------------------------------------------------#
### Histograms for the parameters ###

using DelimitedFiles

# Load the parameters
prms_int = readdlm("WGS_particle_reaction/csv_files/samples_LHS_init.csv", ',', Float64, '\n')
prms_beyond = readdlm("WGS_particle_reaction/csv_files/samples_LHS_beyond.csv", ',', Float64, '\n')
prms_random = readdlm("WGS_particle_reaction/csv_files/samples_LHS_random.csv", ',', Float64, '\n')

using Plots

## Init ##

histogram(prms_int[:, 1], label = "T", dpi = 600)
title!("Histogram of the T for the first approach")
xlabel!("Temperature [K]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/T_init.png")

histogram(prms_int[:, 2], label = "C_1", dpi = 600)
title!("Histogram of the C_1 for the first approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C1_init.png")

histogram(prms_int[:, 3], label = "C_2", dpi = 600)
title!("Histogram of the C_2 for the first approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C2_init.png")

histogram(prms_int[:, 4], label = "C_3", dpi = 600)
title!("Histogram of the C_3 for the first approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C3_init.png")

histogram(prms_int[:, 5], label = "C_4", dpi = 600)
title!("Histogram of the C_4 for the first approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C4_init.png")

histogram(prms_int[:, 6], label = "C_5", dpi = 600)
title!("Histogram of the C_5 for the first approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C5_init.png")

## Beyond ##

histogram(prms_beyond[:, 1], label = "T", dpi = 600)
title!("Histogram of the T for the second approach")
xlabel!("Temperature [K]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/T_beyond.png")

histogram(prms_beyond[:, 2], label = "C_1", dpi = 600)
title!("Histogram of the C_1 for the second approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C1_beyond.png")

histogram(prms_beyond[:, 3], label = "C_2", dpi = 600)
title!("Histogram of the C_2 for the second approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C2_beyond.png")

histogram(prms_beyond[:, 4], label = "C_3", dpi = 600)
title!("Histogram of the C_3 for the second approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C3_beyond.png")

histogram(prms_beyond[:, 5], label = "C_4", dpi = 600)
title!("Histogram of the C_4 for the second approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C4_beyond.png")

histogram(prms_beyond[:, 6], label = "C_5", dpi = 600)
title!("Histogram of the C_5 for the second approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C5_beyond.png")

## Random ##

histogram(prms_random[:, 1], label = "T", dpi = 600)
title!("Histogram of the T for the third approach")
xlabel!("Temperature [K]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/T_random.png")

histogram(prms_random[:, 2], label = "C_1", dpi = 600)
title!("Histogram of the C_1 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C1_random.png")

histogram(prms_random[:, 3], label = "C_2", dpi = 600)
title!("Histogram of the C_2 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C2_random.png")

histogram(prms_random[:, 4], label = "C_3", dpi = 600)
title!("Histogram of the C_3 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C3_random.png")

histogram(prms_random[:, 5], label = "C_4", dpi = 600)
title!("Histogram of the C_4 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C4_random.png")

histogram(prms_random[:, 6], label = "C_5", dpi = 600)
title!("Histogram of the C_5 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/C5_random.png")

histogram(prms_random[:, 7], label = "C_c_1", dpi = 600)
title!("Histogram of the C_c_1 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/Cc1_random.png")

histogram(prms_random[:, 8], label = "C_c_2", dpi = 600)
title!("Histogram of the C_c_2 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/Cc2_random.png")

histogram(prms_random[:, 9], label = "C_c_3", dpi = 600)
title!("Histogram of the C_c_3 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/Cc3_random.png")

histogram(prms_random[:, 10], label = "C_c_4", dpi = 600)
title!("Histogram of the C_c_4 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/Cc4_random.png")

histogram(prms_random[:, 11], label = "C_c_5", dpi = 600)
title!("Histogram of the C_c_5 for the third approach")
xlabel!("Concentration [mol/m^3]")
ylabel!("Frequency")
savefig("WGS_particle_reaction/appx_figures/Cc5_random.png")

# ------------------------------------------------------------------------------------------------------------------------------------------- #
using DelimitedFiles
using Statistics

data_m2 = readdlm("WGs_particle_reaction/csv_files/all_data_m2.csv", ',', Float64, '\n')[:,12:16]

m2_c1 = abs.(data_m2[:, 1])
m2_c2 = abs.(data_m2[:, 2])
m2_c3 = abs.(data_m2[:, 3])
m2_c4 = abs.(data_m2[:, 4])
m2_c5 = abs.(data_m2[:, 5])

min_m2 = [minimum(m2_c1) minimum(m2_c2) minimum(m2_c3) minimum(m2_c4) minimum(m2_c5)]
max_m2 = [maximum(m2_c1) maximum(m2_c2) maximum(m2_c3) maximum(m2_c4) maximum(m2_c5)]

data_m2_test = readdlm("WGs_particle_reaction/csv_files/test_data_m2.csv", ',', Float64, '\n')[:,12:16]

m2_c1_test = abs.(data_m2_test[:, 1])
m2_c2_test = abs.(data_m2_test[:, 2])
m2_c3_test = abs.(data_m2_test[:, 3])
m2_c4_test = abs.(data_m2_test[:, 4])
m2_c5_test = abs.(data_m2_test[:, 5])

min_m2_test = [minimum(m2_c1_test) minimum(m2_c2_test) minimum(m2_c3_test) minimum(m2_c4_test) minimum(m2_c5_test)]
max_m2_test = [maximum(m2_c1_test) maximum(m2_c2_test) maximum(m2_c3_test) maximum(m2_c4_test) maximum(m2_c5_test)]
# ------------------------------------------------------------------------------------------------------------------------------------------- #
using DelimitedFiles
using Plots

k1_t = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_sol_t_1e-6.csv", ',', Float64, '\n')*3600

k1_dcc1 = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_dCC1dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k1_dcc2 = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_dCC2dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k1_dcc3 = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_dCC3dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k1_dcc4 = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_dCC4dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k1_dcc5 = readdlm("WGS_particle_reaction/ml_data_all/k1/k1_dCC5dt_1e-6.csv", ',', Float64, '\n')[:, 1]

plot(k1_t, k1_dcc1, label = "dC_c_1/dt", dpi = 600)
plot!(k1_t, k1_dcc2, label = "dC_c_2/dt")
plot!(k1_t, k1_dcc3, label = "dC_c_3/dt")
plot!(k1_t, k1_dcc4, label = "dC_c_4/dt")
plot!(k1_t, k1_dcc5, label = "dC_c_5/dt")
title!("Derivative of concentration profiles of C_c_i for k1", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Derivative of concentration [mol/m^3 h]")
savefig("WGS_particle_reaction/appx_figures/k1_rate_of_change.png")

k101_dcc1 = readdlm("WGS_particle_reaction/ml_data_all/k101/k101_dCC1dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k101_dcc2 = readdlm("WGS_particle_reaction/ml_data_all/k101/k101_dCC2dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k101_dcc3 = readdlm("WGS_particle_reaction/ml_data_all/k101/k101_dCC3dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k101_dcc4 = readdlm("WGS_particle_reaction/ml_data_all/k101/k101_dCC4dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k101_dcc5 = readdlm("WGS_particle_reaction/ml_data_all/k101/k101_dCC5dt_1e-6.csv", ',', Float64, '\n')[:, 1]

plot(k1_t, k101_dcc1, label = "dC_c_1/dt", dpi = 600)
plot!(k1_t, k101_dcc2, label = "dC_c_2/dt")
plot!(k1_t, k101_dcc3, label = "dC_c_3/dt")
plot!(k1_t, k101_dcc4, label = "dC_c_4/dt")
plot!(k1_t, k101_dcc5, label = "dC_c_5/dt")
title!("Derivative of concentration profiles of C_c_i for k101", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Derivative of concentration [mol/m^3 h]")
savefig("WGS_particle_reaction/appx_figures/k101_rate_of_change.png")

k10101_dcc1 = readdlm("WGS_particle_reaction/ml_data_all/k10101/k10101_dCC1dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k10101_dcc2 = readdlm("WGS_particle_reaction/ml_data_all/k10101/k10101_dCC2dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k10101_dcc3 = readdlm("WGS_particle_reaction/ml_data_all/k10101/k10101_dCC3dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k10101_dcc4 = readdlm("WGS_particle_reaction/ml_data_all/k10101/k10101_dCC4dt_1e-6.csv", ',', Float64, '\n')[:, 1]
k10101_dcc5 = readdlm("WGS_particle_reaction/ml_data_all/k10101/k10101_dCC5dt_1e-6.csv", ',', Float64, '\n')[:, 1]

plot(k1_t, k10101_dcc1, label = "dC_c_1/dt", dpi = 600)
plot!(k1_t, k10101_dcc2, label = "dC_c_2/dt")
plot!(k1_t, k10101_dcc3, label = "dC_c_3/dt")
plot!(k1_t, k10101_dcc4, label = "dC_c_4/dt")
plot!(k1_t, k10101_dcc5, label = "dC_c_5/dt")
title!("Derivative of concentration profiles of C_c_i for k10101", titlefontsize=10, titlefontcolor=:black)
xlabel!("Time [s]")
ylabel!("Derivative of concentration [mol/m^3 h]")
savefig("WGS_particle_reaction/appx_figures/k10101_rate_of_change.png")
