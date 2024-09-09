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
t0_disc = time()
prob = discretize(WGS_pde, discretization)
t1_disc = time() - t0_disc
print("\nDiscretization time: ", t1_disc, "\n")

# Solving ODE
t0_sol = time()
sol = solve(prob, FBDF(), abstol = 1e-6, reltol = 1e-6)
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

# Sampling ---------------------------------------------------------------------------------------------------------------------------------- #

using Surrogates
using Plots

C_total_max = P_val / (R * 1e-3 * 393)

#           T  y_CO y_CO2  y_H2  y_H2O  y_N2 
p_upper = [593.0, C_total_max, C_total_max, C_total_max, C_total_max, C_total_max]
p_lower = [393.0, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]
nsample = 10000

samples_LHS_og = sample(nsample, p_lower, p_upper, LatinHypercubeSample())
samples_LHS = hcat(map(x -> collect(x), samples_LHS_og)...)
samples_LHS = samples_LHS'

histogram(samples_LHS[:, 1], label = "T")
p1 = histogram(samples_LHS[:, 2], label = "C_CO")
histogram!(p1, samples_LHS[:, 3], label = "C_CO2")
histogram!(p1, samples_LHS[:, 4], label = "C_H2")
histogram!(p1, samples_LHS[:, 5], label = "C_H2O")
histogram!(p1, samples_LHS[:, 6], label = "C_N2")

using DelimitedFiles
write_to_csv("samples_LHS_beyond.csv", samples_LHS, "WGS_particle_reaction")

plot(samples_LHS[:, 1], samples_LHS[:, 2], seriestype = :scatter, label = "LHS samples T vs C_CO")
plot(samples_LHS[:, 1], samples_LHS[:, 3], seriestype = :scatter, label = "LHS samples T vs C_CO2")
plot(samples_LHS[:, 1], samples_LHS[:, 4], seriestype = :scatter, label = "LHS samples T vs C_H2")
plot(samples_LHS[:, 1], samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples T vs C_H2O")
plot(samples_LHS[:, 1], samples_LHS[:, 6], seriestype = :scatter, label = "LHS samples T vs C_N2")
plot(samples_LHS[:, 2], samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples C_CO vs C_H2O")
plot(samples_LHS[:, 1], samples_LHS[:, 2]./samples_LHS[:, 5], seriestype = :scatter, label = "LHS samples T vs syngas ratio")

# ------------------------------------------------------------------------------------------------------------------------------------------- #

using DelimitedFiles

p_sampled = readdlm("WGS_particle_reaction/samples_LHS_beyond.csv", ',', Float64, '\n')

new_cci_init = readdlm("WGS_particle_reaction/cci_t0_beyond.csv", ',', Float64, '\n')

time_loop = []
error_loop = Any[("k", "p")]

# --------------------------------------------------------------#
# to make sure that k and p are not reset after inner loop

# put final k number here (from last run)
k_num_0 = 100
k_num = 100

# put final p number here (p_num = k_num - 100)
p_num = k_num - 100

# j should increase by 1 every 100 k (with current setup, change if needed)
j_num = floor(Int, k_num / 100)
# --------------------------------------------------------------#

t0_total = time()

for j in (j_num):(j_num + 2)

    C_c_i_init_new = new_cci_init[j, :]

    for k in 1:100
        t0_it = time()

        k_num = k_num + 1
        p_num = p_num + 1

        folder = "WGS_particle_reaction/ml_data_beyond/k$(k_num)"

        try
            newprms_scal = [T => p_sampled[p_num, 1], P => P_val, F => F_val]
            newprms_vec_C_i = [C_i[i] => p_sampled[p_num, i+1] for i in 1:5]
            newprms_vec_C_c_i_init = [C_c_i_init[i] => C_c_i_init_new[i] for i in 1:5]
            newprms = [newprms_scal...; newprms_vec_C_i...; newprms_vec_C_c_i_init...]

            # Domain (time is in [h])
            domains_new = [t ∈ Interval(0.0, 5e-4),
                r ∈ Interval(0.0, rad_cat)]

            @named WGS_pde_new = PDESystem(eqs, bcs, domains_new, [t, r], vars, newprms)
            discretization_new = MOLFiniteDifference([r => dr], t, order=order)

            t0_disc_new = time()
            newprob = discretize(WGS_pde_new, discretization_new)
            t1_disc_new = time() - t0_disc_new
            print("\nDiscretization time k = $(k_num): ", t1_disc_new, "\n")

            t0_sol_new = time()
            newsol = solve(newprob, FBDF(), abstol = 1e-6, reltol = 1e-6)
            t1_sol_new = time() - t0_sol_new
            print("\nSolution time k = $(k_num): ", t1_sol_new, "\n")

            write_to_csv("k$(k_num)_sol_cc1_og.csv", newsol[C_c_1(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc2_og.csv", newsol[C_c_2(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc3_og.csv", newsol[C_c_3(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc4_og.csv", newsol[C_c_4(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc5_og.csv", newsol[C_c_5(t, r)], folder)
            write_to_csv("k$(k_num)_sol_t_og.csv", newsol.t, folder)

            dccidt, _ = finite_diff_sol(newsol)

            write_to_csv("k$(k_num)_dCC1dt_og.csv", dccidt[1], folder)
            write_to_csv("k$(k_num)_dCC2dt_og.csv", dccidt[2], folder)
            write_to_csv("k$(k_num)_dCC3dt_og.csv", dccidt[3], folder)
            write_to_csv("k$(k_num)_dCC4dt_og.csv", dccidt[4], folder)
            write_to_csv("k$(k_num)_dCC5dt_og.csv", dccidt[5], folder)

            plot(newsol.t, newsol[C_c_1(t, r)], label = "CO")
            plot!(newsol.t, newsol[C_c_2(t, r)], label = "CO2")
            plot!(newsol.t, newsol[C_c_3(t, r)], label = "H2")
            plot!(newsol.t, newsol[C_c_4(t, r)], label = "H2O")
            plot!(newsol.t, newsol[C_c_5(t, r)], label = "N2")
            savefig("WGS_particle_reaction/ml_data_beyond/k$(k_num)/k$(k_num)_sol_og.png")

            plot(newsol.t, dccidt[1], label = "dCC_CO_dt")
            plot!(newsol.t, dccidt[2], label = "dCC_CO2_dt")
            plot!(newsol.t, dccidt[3], label = "dCC_H2_dt")
            plot!(newsol.t, dccidt[4], label = "dCC_H2O_dt")
            plot!(newsol.t, dccidt[5], label = "dCC_N2_dt")
            savefig("WGS_particle_reaction/ml_data_beyond/k$(k_num)/k$(k_num)_dCCdt_og.png")

            # newsol with constant timestep
            newsol_const = solve(newprob, FBDF(), saveat = 1e-6 ,abstol = 1e-6, reltol = 1e-6)

            write_to_csv("k$(k_num)_sol_cc1_1e-6.csv", newsol_const[C_c_1(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc2_1e-6.csv", newsol_const[C_c_2(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc3_1e-6.csv", newsol_const[C_c_3(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc4_1e-6.csv", newsol_const[C_c_4(t, r)], folder)
            write_to_csv("k$(k_num)_sol_cc5_1e-6.csv", newsol_const[C_c_5(t, r)], folder)
            write_to_csv("k$(k_num)_sol_t_1e-6.csv", newsol_const.t, folder)

            dccidt_const, _ = finite_diff_sol(newsol_const)

            write_to_csv("k$(k_num)_dCC1dt_1e-6.csv", dccidt_const[1], folder)
            write_to_csv("k$(k_num)_dCC2dt_1e-6.csv", dccidt_const[2], folder)
            write_to_csv("k$(k_num)_dCC3dt_1e-6.csv", dccidt_const[3], folder)
            write_to_csv("k$(k_num)_dCC4dt_1e-6.csv", dccidt_const[4], folder)
            write_to_csv("k$(k_num)_dCC5dt_1e-6.csv", dccidt_const[5], folder)

            plot(newsol_const.t, newsol_const[C_c_1(t, r)], label = "CO")
            plot!(newsol_const.t, newsol_const[C_c_2(t, r)], label = "CO2")
            plot!(newsol_const.t, newsol_const[C_c_3(t, r)], label = "H2")
            plot!(newsol_const.t, newsol_const[C_c_4(t, r)], label = "H2O")
            plot!(newsol_const.t, newsol_const[C_c_5(t, r)], label = "N2")
            savefig("WGS_particle_reaction/ml_data_beyond/k$(k_num)/k$(k_num)_sol_1e-6.png")

            plot(newsol_const.t, dccidt_const[1], label = "dCC_CO_dt")
            plot!(newsol_const.t, dccidt_const[2], label = "dCC_CO2_dt")
            plot!(newsol_const.t, dccidt_const[3], label = "dCC_H2_dt")
            plot!(newsol_const.t, dccidt_const[4], label = "dCC_H2O_dt")
            plot!(newsol_const.t, dccidt_const[5], label = "dCC_N2_dt")
            savefig("WGS_particle_reaction/ml_data_beyond/k$(k_num)/k$(k_num)_dCCdt_1e-6.png")
        catch
            print("\nError with k = $(k_num), p = $(p_num) \n")
            push!(error_loop, (k_num, p_num))
        end
        
        t1_it = time() - t0_it
        push!(time_loop, t1_it)

        print("\nDone with k = $(k_num), p = $(p_num) \n")
    end
end

t1_total = time() - t0_total
push!(time_loop, t1_total)
print("\nTotal time: ", t1_total, "\n")

k_num_0 = k_num_0 + 1
write_to_csv("time_k$(k_num_0)-k$(k_num).csv", time_loop, "WGS_particle_reaction/ml_data_beyond")

# Document failed runs ------------------------------------------------------------------------------- #
write_to_csv("errors_k$(k_num_0)-k$(k_num).csv", error_loop, "WGS_particle_reaction/ml_data_beyond")
# ---------------------------------------------------------------------------------------------------- #

print("final p_num: ", p_num, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------- #

