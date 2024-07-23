#### Test functions for the WGS model (mortly outdated) ####

using Pkg
Pkg.activate("WGS")
# include("functions_WGS.jl")

### Test functions ###

## Test μ_mix_func ##
y = [0.5, 0.3, 0.2]
μ_i= [11.1e-6, 9.0e-6, 8.0e-6]
M = [16.04, 30.07, 44.10]

μ_mix = μ_mix_func(y, μ_i, M)

## Test D_eff_ij_func and D_i_m_func ##
theta = 0.8
tau = 1.0
y = [0.2, 0.2, 0.2, 0.2, 0.2]
# D_ij_mat = D_ij_matrix_func(300, 1)

#D_eff_ij = zeros(size(D_ij_mat))

D_eff_ij = D_eff_ij_func(D_ij_mat, theta, tau)

# D_i_m = D_i_m_func(y, D_eff_ij)


## Test C_p_func ##
y = [0.5, 0.3, 0.2]
C_p_i = [28.1, 33.2, 42.0] 

C_p_mix = C_p_func(y, C_p_i)

## Test A_ij_func, λ_dash_func ##
y_list = [0.23, 0.77]
T_example = 373.15
k_list = [0.02504, 0.01587]
μ_list = [1.161e-5, 1.361e-5]
M_list = [46.07, 50.49]
T_boil_list = [248.3, 248.9]

A_11 = A_ij_func(1, 1, μ_list, M_list, T_example, T_boil_list, 1.0)
A_12 = A_ij_func(1, 2, μ_list, M_list, T_example, T_boil_list, 1.0)
A_21 = A_ij_func(2, 1, μ_list, M_list, T_example, T_boil_list, 1.0)
A_22 = A_ij_func(2, 2, μ_list, M_list, T_example, T_boil_list, 1.0)

λ_dash = λ_dash_func(y_list, k_list, μ_list, M_list, T_example, T_boil_list, 1.0)

## Test λ_func ## 
M = 44.01
P_c = 7.383
V_c = 0.0940
Z_c = P_c * V_c / (0.008314 * T_c)
# Z_c1 = 0.274
k_dash = 0.0220
V = 0.22809
ρ = V_c/V
# ρ1 = 0.4121

λ_test = λ_func(y, T, P, R, M, k_dash)

## Test syntax/use of symbolic functions with array as ouput ##
using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets

function test(u)
    return u ./ 5
end


@register_symbolic test(u)

n_comp = 2
# @register_array_symbolic test(u) begin
#     size=(n_comp,)
#     eltype=eltype(Vector{Num})
# end

# Parameters, variables, and derivatives

@parameters t, x, p[1:n_comp], q[1:n_comp] , f
@variables T[1:n_comp], u(..)[1:n_comp]
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
params = Symbolics.scalarize(reduce(vcat,[p .=> [1.5, 2.0], q .=> [1.2, 1.8], f => 2]))
# 1D PDE and boundary conditions

#eqs = [Dt(u(t, x)[i]) ~ p[i] * Dxx(u(t, x)[i]) for i in 1:n_comp]
#eqs1 = [T[i] ~ test(u(t, x))[i] for i in 1:n_comp]
eqs = [Dt(u(t, x)[i]) ~ f * p[i] * test(u(t, x)[i]) * Dxx(u(t, x)[i]) for i in 1:n_comp]
#eqs = [eqs1; eqs2]
bcs = [[u(0, x)[i] ~ q[i] * cos(x),
        u(t, 0)[i] ~ sin(t),
        u(t, 1)[i] ~ exp(-t) * cos(1),
        Dx(u(t,0)[i]) ~ 0.0] for i in 1:n_comp]
bcs_collected = reduce(vcat, bcs)

# Space and time domains
domains = [t ∈ Interval(0.0, 1.0),
           x ∈ Interval(0.0, 1.0)]

# PDE system

@named pdesys = PDESystem(eqs, bcs_collected, domains, [t, x], [u(t, x)[1:n_comp]], params)


# Method of lines discretization
dx = 0.1
order = 2
discretization = MOLFiniteDifference([x => dx], t; approx_order = order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization) #error occurs here

# Solve ODE problem
sol = solve(prob, Tsit5(), saveat=0.2)

# for derivatives of functions:
#= @parameters r, t, z
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
@variables P(t, r, z), T(z, t, r), C(r, z), (C_i(z, t))[1:5]

teste = Dr.(D_ij_matrix_func(T, P, matrix_D_ij))
derivatives = expand_derivatives.(teste)
## Register symbolic functions ###
derivatives[1, 2] =#

# y_i = C_i_val / sum(C_i_val)
# n = F_0 * y_i
# V = π * (D_rct_val/2)^2 * L_val
# C = n / V
# C./C_i_val


# f1 = 1e-5
# n1 = f1 * y_i
# V1 = π * ((12.7e-3)/2)^2 * (4.8e-3)
# c1 = n1 / V1

C_c_i_init = [0.0, 0.0, 0.0, 0.0, 16.44604367824988]