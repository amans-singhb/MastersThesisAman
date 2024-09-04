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

## Constant parameters

# Dimensions reactor
D_rct = 12.7e-3 # [m] 
L_rct = 4.8e-3 # [m]
V_rct = π * (D_rct/2)^2 * L_rct # [m^3]

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
f_temp = F_ind * volume_ratio
# ----------------------------------------------------------------------------- #
F_val = f_temp * (P_val/ P_ind) * (T_ind / T_val) # [mol/h] [param3]

# Temp vals C ----------------------------------------------------------------- #
C_total = P_val / (R * 1e-3 * T_val) # [mol/m^3]
C_i_val_temp = C_total * y_0 # [mol/m^3]
C_i_val_temp[4]/C_i_val_temp[3] # 2.3
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
domains = [t ∈ Interval(0.0, 3e-3),
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
t0_disc = time()
prob = discretize(WGS_pde, discretization)
t1_disc = time() - t0_disc
print("\nDiscretization time: ", t1_disc, "\n")

# Solving ODE
t0_sol = time()
sol = solve(prob, FBDF(), saveat = 1e-7 ,abstol = 1e-6, reltol = 1e-6)
t1_sol = time() - t0_sol
print("\nSolution time: ", t1_sol, "\n")

using Plots

sols1 = sol[C_c_1(t, r)][:, 21]
sols2 = sol[C_c_2(t, r)][:, 21]
sols3 = sol[C_c_3(t, r)][:, 21]
sols4 = sol[C_c_4(t, r)][:, 21]
sols5 = sol[C_c_5(t, r)][:, 21]
sol_t = sol.t * 3600

plot(sol_t, sols1, label = "CO")
plot!(sol_t, sols2, label = "CO2")
plot!(sol_t, sols3, label = "H2")
plot!(sol_t, sols4, label = "H2O")
plot!(sol_t, sols5, label = "N2")
# savefig("WGS_particle_reaction/sol_postdiff2.png")

# New parameters ---------------------------------------------------------------------------------------------------------------------------- #
C_i_old = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]

# postdiff
C_c_i_init_new = C_i_old

# # prediff
# zero_val = 1e-15 # isapprox(0, 1e-324) = true
# C_c_i_init_new = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]

newprms_scal = [T => 483, P => 2, F => 2e-5]
newprms_vec_C_i = [C_i[i] => C_i_old[i] for i in 1:5]
newprms_vec_C_c_i_init = [C_c_i_init[i] => C_c_i_init_new[i] for i in 1:5]
newprms = [newprms_scal...; newprms_vec_C_i...; newprms_vec_C_c_i_init...]

# # Domain postdiff (time is in [h])
# domains_new = [t ∈ Interval(0.0, 2e-5),
#     r ∈ Interval(0.0, rad_cat)]

# Domain prediff (time is in [h])
domains_new = [t ∈ Interval(0.0, 3e-3),
    r ∈ Interval(0.0, rad_cat)]

@named WGS_pde_new = PDESystem(eqs, bcs, domains_new, [t, r], vars, newprms)
discretization_new = MOLFiniteDifference([r => dr], t, order=order)

t0_disc_new = time()
newprob = discretize(WGS_pde_new, discretization_new)
t1_disc_new = time() - t0_disc_new
print("\nDiscretization time: ", t1_disc_new, "\n")

t0_sol_new = time()
newsol = solve(newprob, KenCarp47(), saveat = 1e-7 , abstol = 1e-6, reltol = 1e-6)
t1_sol_new = time() - t0_sol_new
print("\nSolution time: ", t1_sol_new, "\n")

newsols1 = newsol[C_c_1(t, r)][:, 21]
newsols2 = newsol[C_c_2(t, r)][:, 21]
newsols3 = newsol[C_c_3(t, r)][:, 21]
newsols4 = newsol[C_c_4(t, r)][:, 21]
newsols5 = newsol[C_c_5(t, r)][:, 21]
newsol_t = newsol.t * 3600

plot(newsol_t, newsols1, label = "CO")
plot!(newsol_t, newsols2, label = "CO2")
plot!(newsol_t, newsols3, label = "H2")
plot!(newsol_t, newsols4, label = "H2O")
plot!(newsol_t, newsols5, label = "N2")
# savefig("WGS_particle_reaction/sol_prediff2.png")

# Numerical differentiation ---------------------------------------------------------------------------------------------------------------- #

dCCidt, error_m  = finite_diff_sol(sol)

plot(sol.t, error_m[1], label = "error_CO")

r_idx = 1
plot(sol.t, dCCidt[1][:, r_idx], label = "dCC_CO_dt")
plot!(sol.t, dCCidt[2][:, r_idx], label = "dCC_CO2_dt")
plot!(sol.t, dCCidt[3][:, r_idx], label = "dCC_H2_dt")
plot!(sol.t, dCCidt[4][:, r_idx], label = "dCC_H2O_dt")
plot!(sol.t, dCCidt[5][:, r_idx], label = "dCC_N2_dt")

plot(sol.t, dCCidt[1][:, 1], label = "dCC_CO_dt1")
plot!(sol.t, dCCidt[1][:, 21], label = "dCC_CO_dt21")

plot(sol.t, dCCidt[2][:, 1], label = "dCC_CO2_dt1")
plot!(sol.t, dCCidt[2][:, 21], label = "dCC_CO2_dt21")

plot(sol.t, dCCidt[3][:, 1], label = "dCC_H2_dt1")
plot!(sol.t, dCCidt[3][:, 21], label = "dCC_H2_dt21")

plot(sol.t, dCCidt[4][:, 1], label = "dCC_H2O_dt1")
plot!(sol.t, dCCidt[4][:, 21], label = "dCC_H2O_dt21")

plot(sol.t, dCCidt[5][:, 1], label = "dCC_N2_dt1")
plot!(sol.t, dCCidt[5][:, 21], label = "dCC_N2_dt21")

# Numerical integration -------------------------------------------------------------------------------------------------------------------- #

dCC1dt = dCCidt[1]
sol_cc1 = sol[C_c_1(t, r)]

dCC2dt = dCCidt[2]
sol_cc2 = sol[C_c_2(t, r)]

dCC3dt = dCCidt[3]
sol_cc3 = sol[C_c_3(t, r)]

dCC4dt = dCCidt[4]
sol_cc4 = sol[C_c_4(t, r)]

dCC5dt = dCCidt[5]
sol_cc5 = sol[C_c_5(t, r)]

# Implicit Euler
dt = 1e-7
time_len = 0.0:dt:3e-3

CC1 = zeros(size(sol_cc1))
CC2 = zeros(size(sol_cc2))
CC3 = zeros(size(sol_cc3))
CC4 = zeros(size(sol_cc4))
CC5 = zeros(size(sol_cc5))

CC1[1,:] = [C_c_i_init_val[1] for i in 1:21]
CC2[1,:] = [C_c_i_init_val[2] for i in 1:21]
CC3[1,:] = [C_c_i_init_val[3] for i in 1:21]
CC4[1,:] = [C_c_i_init_val[4] for i in 1:21]
CC5[1,:] = [C_c_i_init_val[5] for i in 1:21]


for i in 2:length(time_len)
    CC1[i, :] = CC1[i-1, :] + dt .* dCC1dt[i, :]
    CC2[i, :] = CC2[i-1, :] + dt .* dCC2dt[i, :]
    CC3[i, :] = CC3[i-1, :] + dt .* dCC3dt[i, :]
    CC4[i, :] = CC4[i-1, :] + dt .* dCC4dt[i, :]
    CC5[i, :] = CC5[i-1, :] + dt .* dCC5dt[i, :]
end

r_index = 20

plot(time_len, CC1[:, r_index], label = "CO (implicit Euler), r = $r_index")
plot!(sol.t, sol_cc1[:, r_index], label = "CO (true), r = $r_index")

plot(time_len, CC2[:, r_index], label = "CO2 (implicit Euler), r = $r_index")
plot!(sol.t, sol_cc2[:, r_index], label = "CO2 (true), r = $r_index")

plot(time_len, CC3[:, r_index], label = "H2 (implicit Euler), r = $r_index")
plot!(sol.t, sol_cc3[:, r_index], label = "H2 (true), r = $r_index")

plot(time_len, CC4[:, r_index], label = "H2O (implicit Euler), r = $r_index")
plot!(sol.t, sol_cc4[:, r_index], label = "H2O (true), r = $r_index")

plot(time_len, CC5[:, r_index], label = "N2 (implicit Euler), r = $r_index")
plot!(sol.t, sol_cc5[:, r_index], label = "N2 (true), r = $r_index")

# Error calculation
error_cc1 = abs.(CC1 - sol_cc1)
error_cc2 = abs.(CC2 - sol_cc2)
error_cc3 = abs.(CC3 - sol_cc3)
error_cc4 = abs.(CC4 - sol_cc4)
error_cc5 = abs.(CC5 - sol_cc5)

plot(time_len, error_cc1[:, r_index], label = "error CO (implicit Euler), r = $r_index")
plot!(time_len, error_cc2[:, r_index], label = "error CO2 (implicit Euler), r = $r_index")
plot!(time_len, error_cc3[:, r_index], label = "error H2 (implicit Euler), r = $r_index")
plot!(time_len, error_cc4[:, r_index], label = "error H2O (implicit Euler), r = $r_index")
plot!(time_len, error_cc5[:, r_index], label = "error N2 (implicit Euler), r = $r_index")

plot_cc1_err =  plot(time_len, error_cc1[:, 1], label = "error CO (implicit Euler), r = 1")
for i in 2:21
    plot!(plot_cc1_err, time_len, error_cc1[:, i], label = "error CO (implicit Euler), r = $i")
end
display(plot_cc1_err)

plot_cc2_err =  plot(time_len, error_cc2[:, 1], label = "error CO2 (implicit Euler), r = 1")
for i in 2:21
    plot!(plot_cc2_err, time_len, error_cc2[:, i], label = "error CO2 (implicit Euler), r = $i")
end
display(plot_cc2_err)

plot_cc3_err =  plot(time_len, error_cc3[:, 1], label = "error H2 (implicit Euler), r = 1")
for i in 2:21
    plot!(plot_cc3_err, time_len, error_cc3[:, i], label = "error H2 (implicit Euler), r = $i")
end
display(plot_cc3_err)

plot_cc4_err =  plot(time_len, error_cc4[:, 1], label = "error H2O (implicit Euler), r = 1")
for i in 2:21
    plot!(plot_cc4_err, time_len, error_cc4[:, i], label = "error H2O (implicit Euler), r = $i")
end
display(plot_cc4_err)

plot_cc5_err =  plot(time_len, error_cc5[:, 1], label = "error N2 (implicit Euler), r = 1")
for i in 2:21
    plot!(plot_cc5_err, time_len, error_cc5[:, i], label = "error N2 (implicit Euler), r = $i")
end
display(plot_cc5_err)
#------------------------------------------------------------------------------------------------------------------------------------------- #

using Surrogates
#           T  y_CO y_CO2  y_H2  y_H2O  y_N2 
p_upper = [593.0, 1.0, 1.0, 1.0, 1.0, 1.0]
p_lower = [393.0, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]
nsample = 2000

samples_LHS_og = sample(nsample, p_lower, p_upper, LatinHypercubeSample())
samples_LHS = hcat(map(x -> collect(x), samples_LHS_og)...)
samples_LHS = samples_LHS'

# Ensuring the molcefractions sum up to 1
samples_LHS[:, 2:end] = samples_LHS[:, 2:end] ./ sum(samples_LHS[:, 2:end], dims = 2)
for i in 1:nsample
    if !isapprox(1, sum(samples_LHS[i, 2:end]))
        println("Error in sample $i")
    end
end
samples_LHS

# Converting molcefractions to concentrations based on new temperature
for i in 1:nsample
    C_tot_in  = P_val / (R * 1e-3 * samples_LHS[i, 1])
    samples_LHS[i, 2:end] = samples_LHS[i, 2:end] .* C_tot_in
end

samples_LHS
for i in 1:nsample
    C_tot_in  = P_val / (R * 1e-3 * samples_LHS[i, 1])
    if C_tot_in != sum(samples_LHS[i, 2:end])
        println("1. Error in sample $i, difference= ", C_tot_in - sum(samples_LHS[i, 2:end]))
    end
end

using Plots

plot(samples_LHS[:, 1], samples_LHS[:, 2], seriestype = :scatter, label = "LHS samples T vs C_CO")
plot(samples_LHS[:, 1], samples_LHS[:, 3], seriestype = :scatter, label = "LHS samples T vs C_CO2")
plot(samples_LHS[:, 1], samples_LHS[:, 4], seriestype = :scatter, label = "LHS samples T vs C_H2")
plot(samples_LHS[:, 1], samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples T vs C_H2O")
plot(samples_LHS[:, 1], samples_LHS[:, 6], seriestype = :scatter, label = "LHS samples T vs C_N2")
plot(samples_LHS[:, 2], samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples C_CO vs C_H2O")
plot(samples_LHS[:, 1], samples_LHS[:, 2]./samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples T vs syngas ratio")

