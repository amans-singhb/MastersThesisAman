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

## Parameters

# Inlet values
F_val= 1e-5 # [mol/h] [param3]
T_val = 503.0 # [K] [param1]
P_val = 1.3 # [atm] [param2]

# Temp vals ------------------------------------------------------------------- #
y_0 = [0.208917; 0.0910445; 0.204129; 0.480558; 0.0153515;] # from paper
C_i_val_temp = (F_val * y_0) / V_rct # [mol/m^3]
# ----------------------------------------------------------------------------- #

C_i_val = [C_i_val_temp[1], C_i_val_temp[2], C_i_val_temp[3], 2.3 * C_i_val_temp[1], C_i_val_temp[5]] # [mol/m^3]

# postdiff:
C_c_i_init_val = C_i_val

# # prediff:
# zero_val = 1e-15 # isapprox(0, 1e-324) = true
# C_c_i_init = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]

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

# Domain postdiff (time is in [h])
domains = [t ∈ Interval(0.0, 2e-5),
    r ∈ Interval(0.0, rad_cat)]

# # Domain prediff (time is in [h])
# domains = [t ∈ Interval(0.0, 3e-3),
#     r ∈ Interval(0.0, rad_cat)]

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
sol = solve(prob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
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
savefig("WGS_particle_reaction/sol_postdiff.png")

# New parameters ---------------------------------------------------------------------------------------------------------------------------- #
C_i_old = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]
zero_val = 1e-15 # isapprox(0, 1e-324) = true
C_c_i_init_new = [zero_val, zero_val, zero_val, zero_val, (C_total-4*zero_val)]

newprms_scal = [T => 483, P => 2, F => 2e-5]
newprms_vec_C_i = [C_i[i] => C_i_old[i] for i in 1:5]
newprms_vec_C_c_i_init = [C_c_i_init[i] => C_c_i_init_new[i] for i in 1:5]
newprms = [newprms_scal...; newprms_vec_C_i...; newprms_vec_C_c_i_init...]

domains_new = [t ∈ Interval(0.0, 3e-3),
    r ∈ Interval(0.0, rad_cat)]

@named WGS_pde_new = PDESystem(eqs, bcs, domains_new, [t, r], vars, newprms)
discretization_new = MOLFiniteDifference([r => dr], t, order=order)

t0_disc_new = time()
newprob = discretize(WGS_pde_new, discretization_new)
t1_disc_new = time() - t0_disc_new
print("\nDiscretization time: ", t1_disc_new, "\n")

t0_sol_new = time()
newsol = solve(newprob, KenCarp47(), abstol = 1e-6, reltol = 1e-6)
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
savefig("WGS_particle_reaction/sol_prediff.png")
# ------------------------------------------------------------------------------------------------------------------------------------------- #