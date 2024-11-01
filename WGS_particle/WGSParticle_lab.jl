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

#------------------------------------
L_rct = 4.8e-3 # [m]
V_rct = π * (D_rct_val/2)^2 * L_rct # [m^3]
L_rct2 = 304.8e-3 # [m]
V_rct2 = π * (D_rct_val/2)^2 * L_rct2 # [m^3]

# Industrial reactor parameter from paper ------------------------------------- #
F_ind = 9199 * 1000 # [mol/h]
D_ind = 4.4 # [m]
L_ind = 14 # [m]

V_ind = π * (D_ind/2)^2 * L_ind # [m^3]
P_ind = 54.28 # [atm]
T_ind = 503.0 # [K]

q_ind = (F_ind * R_atmm3 * 1e-3 * T_ind) / P_ind # [m^3/h]
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
F_0 = f_temp * (P_val/ P_ind) * (T_ind / T_val) # [mol/h] [param3]

# Temp vals C ----------------------------------------------------------------- #
C_total = P_val / (R_atmm3 * 1e-3 * T_val) # [mol/m^3]
C_i_val_temp = C_total * y_0 # [mol/m^3]
C_i_val_temp[4]/C_i_val_temp[1] # 2.3
# ----------------------------------------------------------------------------- #
C_i_val = [C_i_val_temp[1], C_i_val_temp[2], C_i_val_temp[3], 2.3 * C_i_val_temp[1], C_i_val_temp[5]] # [mol/m^3]

# # postdiff:
# C_c_i_init_val = C_i_val

# prediff:
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
domains = [t ∈ Interval(0.0, 5e-4),
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
prob = discretize(WGS_pde, discretization)

sol = solve(prob, FBDF(), saveat = 1e-6, abstol = 1e-6, reltol = 1e-6)

#--------------------------------------------------------------------------------------------------------------#

using DelimitedFiles
using Statistics
using Plots

sols1 = sol[C_c_1(t, r)][1:280, 1]
sols2 = sol[C_c_2(t, r)][1:280, 1]
sols3 = sol[C_c_3(t, r)][1:280, 1]
sols4 = sol[C_c_4(t, r)][1:280, 1]
sols5 = sol[C_c_5(t, r)][1:280, 1]
solt = sol.t * 3600
sol_t = solt[1:280]

plot(sol_t, sols1, label = "C_c_1", dpi = 1000)
plot!(sol_t, sols2, label = "C_c_2")
plot!(sol_t, sols3, label = "C_c_3")
plot!(sol_t, sols4, label = "C_c_4")
plot!(sol_t, sols5, label = "C_c_5")
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
title!("Concentration profiles of C_c_i from the diffusion model at r = 0.0 m", titlefontsize=10, titlefontcolor=:black)
savefig("WGS_particle/result_diff_model.png")

sol1_sss = sol[C_c_1(t, r)][1:280, :]
sol2_sss = sol[C_c_2(t, r)][1:280, :]
sol3_sss = sol[C_c_3(t, r)][1:280, :]
sol4_sss = sol[C_c_4(t, r)][1:280, :]
sol5_sss = sol[C_c_5(t, r)][1:280, :]
sol_t_sss = sol.t[1:280] * 3600

plot(sol_t_sss, sol1_sss[:, 1:20], label = false, dpi = 1000)
plot!(sol_t_sss, sol1_sss[:, 21], label = "C_c_1")
plot!(sol_t_sss, sol2_sss[:, 1:20], label = false)
plot!(sol_t_sss, sol2_sss[:, 21], label = "C_c_2")
plot!(sol_t_sss, sol3_sss[:, 1:20], label = false)
plot!(sol_t_sss, sol3_sss[:, 21], label = "C_c_3")
plot!(sol_t_sss, sol4_sss[:, 1:20], label = false)
plot!(sol_t_sss, sol4_sss[:, 21], label = "C_c_4")
plot!(sol_t_sss, sol5_sss[:, 1:20], label = false)
plot!(sol_t_sss, sol5_sss[:, 21], label = "C_c_5")
xlabel!("Time [s]")
ylabel!("Concentration [mol/m^3]")
title!("Concentration profiles of C_c_i from the diffusion model at all r", titlefontsize=10, titlefontcolor=:black)
savefig("WGS_particle/result_diff_model_all_r.png")

range_R = range(0.0, rad_cat, length = 21)

cc1 = sol[C_c_1(t, r)][end, :]
cc2 = sol[C_c_2(t, r)][end, :]
cc3 = sol[C_c_3(t, r)][end, :]
cc4 = sol[C_c_4(t, r)][end, :]
cc5 = sol[C_c_5(t, r)][end, :]

mean_c = [mean(cc1), mean(cc), mean(cc3), mean(cc4), mean(cc5)]
diff = mean_c - C_i_val

results = [range_R cc1 cc2 cc3 cc4 cc5]