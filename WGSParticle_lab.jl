##### Catalyst particle balance for WGS reactor model #####

using Pkg
Pkg.activate("WGS")
include("functions_WGSParticle.jl")

### Defining naming conventions and numeration ###

# for species i and j have the following numeration:
# 1 - CO
# 2 - CO2
# 3 - H2
# 4 - H2O
# 5 - N2

#### Particle balance WGS reactor model ####

### Listing some assumptions made for this part of the code ###
# 1. T_c is kept constant, as we are only interested in the evolution of C_c_i
# 2. P_c is kept constant

### PDE system ###

# Parameters
D_rct_val = 12.7e-3 # [m] 
rad_cat = 0.125e-3 # [m]
D_cat_val = rad_cat * 2 # [m] 
d_cat_val = 5904 # [kg/m^3] [param5]

θ_val = 0.55 #[-] [param3]
τ_val = 5 #[-] [param4]

F_0 = 1e-5 # [mol/h] 
T_c_val = 500 # [K] [param1]
P_c_val = 1.3 # [atm] [param6]
R_atmm3 = 8.2057e-2 # [m3 atm/kmol K] [param2]

C_i_val = [0.8054052715722035, 4.495365881590822, 4.411693222036165, 6.4630197133702625, 0.2705595896804266]
C_c_i_init = C_i_val * 0.75

# Mass transfer coefficient
D_i_m_bulk = D_i_m_func(C_i_val, θ_val, τ_val, T_c_val, P_c_val) # [m^2/h]
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
#eqs_ri = [scalarize(r_i .~ r_i_func(C_c_i, d_cat, θ, P_c, T_c))...]
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
domains = [t ∈ Interval(0.0, 0.01),
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
sol = solve(prob, KenCarp47(), saveat = 0.0001, abstol = 1e-6, reltol = 1e-6)
# sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
sols = sol[C_c_1(t, r)]

# Plotting 
time = 0.01
index_sol = Int(time/0.0001)
solution = sols[1:index_sol, :]

using Plots

plot(solution)

using DelimitedFiles
folder = "results_particle_lab"
write_to_csv("C_c_1_lab.csv", sol[C_c_1(t, r)], folder)
write_to_csv("C_c_2_lab.csv", sol[C_c_2(t, r)], folder)
write_to_csv("C_c_3_lab.csv", sol[C_c_3(t, r)], folder)
write_to_csv("C_c_4_lab.csv", sol[C_c_4(t, r)], folder)
write_to_csv("C_c_5_lab.csv", sol[C_c_5(t, r)], folder)

#### Particle balance for diffusion within a sphere ####

### Listing some assumptions made for this part of the code ###
# 1. constant diffusivity
# 2. no Reaction
# 3. constant temperature and pressure

function C_c_i_diff_sphere(M_i, radius, D_i, t, r, N)
    sum_term = 0
    for n in 1:N
        exponent_term = (- n^2 * π^2 * D_i * t) / radius^2
        exp_term = exp(exponent_term)
        sinus_term = (π * n * r) / radius
        sin_term = sin(sinus_term)
        sum_term += n/r * sin_term * exp_term
        print("n = ", n, ": \t" , sum_term,", ", sin_term, ", ", sinus_term, ", ", exp_term, ", ", exponent_term, "\n")
    end
    
    c = (M_i / (2 * radius^2)) * sum_term
    
    return c
end

t_span = range(0, 0.01, length = 101)
r_span = range(0, rad_cat, length = 21)

y_1 = C_c_i_init / sum(C_c_i_init)
n_1 = F_0 * y_1
D_1 = sol[D_1_m]

D_1[2,2], t_span[2], r_span[2]

c_1 = C_c_i_diff_sphere(n_1[1], rad_cat, D_1[2,2], t_span[1], r_span[2], 1)

c_2 = C_c_i_diff_sphere(C_c_i_init[1], rad_cat, D_i_m_bulk[1], t_span[1], r_span[2], 1)
