#### MWE for WGS reactor model ####

### Define functions ###
using Pkg
Pkg.activate("WGS")
# Pkg.add("ModelingToolkit")
# Pkg.add("Symbolics")
# Pkg.add("OrdinaryDiffEq")
# Pkg.add("DomainSets")
# Pkg.add("MethodOfLines"))
using Symbolics
using ModelingToolkit

# Binary gas diffusivities for component pairs {T [K], P [atm], D_ij [cm^2/s]}
function D_ij_func_a(T, P, A, B, C, D, E, F)
    (A * T^B / P) * (log(C / T))^(-2 * D) * exp(-E / T - F / T^2)
end

function D_ij_func_b(P, B)
    B / P
end

function D_ij_func_c(T, P, A, B)
    (A * T + B) / P
end

# Matrix of all binary gas diffusivities for component pairs (SPECIFIC TO THE SYSTEM)
function D_ij_matrix_func(T, P, D_ij_matrix) 
    #p = [eq, i, j, A,      B,      C,  D,  E,   F] (C is set to 1 where it has no value, to avoid log(0) error, D accounts for lack of C in eq = "a" )
    p = ["a" 1 2 3.15e-5 1.57 1 0 113.6 0]
    
    for row in eachrow(p)
        i = row[2]
        j = row[3]

        if row[1] == "a"
            D_ij_matrix[i, j] = D_ij_func_a(T, P, row[4], row[5], row[6], row[7], row[8], row[9])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        elseif row[1] == "b"
            D_ij_matrix[i, j] = D_ij_func_b(P, row[5])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        elseif row[1] == "c"
            D_ij_matrix[i, j] = D_ij_func_c(T, P, row[4], row[5])
            D_ij_matrix[j, i] = D_ij_matrix[i, j]
        else
            print("Invalid equation! Eq: ", row[1])
            return
        end
    end

    return D_ij_matrix
end

# Effective diffusivity of i in j ###[TESTED]###
function D_eff_ij_func(D_ij, θ, τ)
    D_ij .* (θ / τ)
end

# Effective diffusivity of i in the mixture ###[TESTED]###
function D_i_m_func(c, θ, τ, T, P)
    y = [c[1]/(c[1] + c[2]); c[2]/(c[1] + c[2])]
    D_ij_matrix = zeros(Num, 2, 2)
    D_i_m_vec = zeros(Num, 2)

    D_ij = D_ij_matrix_func(T, P, D_ij_matrix)
    D_eff_ij = D_eff_ij_func(D_ij, θ, τ)
    
    for i in eachindex(y)
        i = i[1]
        denominator = 0
        row = D_eff_ij[i, :]

        for j in eachindex(y)
            j = j[1]
            if j != i
                denominator += y[j] / row[j]
            end
        end

        D_i_m_vec[i] = (1-y[i]) / denominator
    end

    D_i_m_vec = D_i_m_vec * 10000 # for conversion

    return D_i_m_vec
end

### PDE system ###

D_cat_val = 0.25e-1 # [m]
rad_cat = 0.5 * D_cat_val # [m]
L_val = 4.8e-3 # [m]

T_val = 450.0 # [K]
T_c_val = 300.0 # [K]
C_i_val = [5.0, 7.0]
drTc_val = 5
θ_val = 0.55 #[-]
τ_val = 5 #[-]
P_c_val = 1.3 # [atm]

C_c_i_init = [0.001, 0.001]

## Parameters ##
@parameters t r T_c θ τ P_c C_i[1:2]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
#@variables C_c_i(..)[1:2] y_c(t, r)[1:2]
@variables C_c_1(..) C_c_2(..)


## Equations and Differential Equations ##
## Differential of function ##

#eqs_y = [y_c[i] ~ C_c_i(t, r)[i] / (C_c_i(t, r)[1] + C_c_i(t, r)[2]) for i in 1:2]

Dr_D_im = Dr.(D_i_m_func([C_c_1(t, r), C_c_2(t, r)], θ, τ, T_c, P_c))
expand_Dr_D_im = expand_derivatives.(Dr_D_im)

eqs1 = [Dt(C_c_1(t, r)) ~ expand_Dr_D_im[1] * Dr(C_c_1(t, r)) + Drr(C_c_1(t, r)),  
Dt(C_c_2(t, r)) ~ expand_Dr_D_im[2] * Dr(C_c_2(t, r)) + Drr(C_c_2(t, r))]
eqs = [eqs1...]

ICS_C_c_i = [C_c_1(0.0, r) ~ C_c_i_init[1], C_c_2(0.0, r) ~ C_c_i_init[2]]

BCS2 = [Dr(C_c_1(t, 0)) ~ 0.0, Dr(C_c_2(t, 0)) ~ 0.0]
BCS3 = [(C_c_1(t, rad_cat) - C_i[1]) ~ 0.0, (C_c_2(t, rad_cat) - C_i[2]) ~ 0.0]


bcs = [ICS_C_c_i...; BCS2...; BCS3...]

using OrdinaryDiffEq, DomainSets, MethodOfLines
using ModelingToolkit: scalarize


# Domain
domains = [t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_c_1(t, r), C_c_2(t, r)]
params_scal = [T_c => T_c_val, θ => θ_val, τ => τ_val, P_c => P_c_val]
params_vec = [C_i[i] => C_i_val[i] for i in 1:2]
params = [params_scal...; params_vec...]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = rad_cat/20
order = 4
discretization = MOLFiniteDifference([r => dr], t, order = order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)
sol = solve(prob, FBDF(), saveat = 0.001, abstol = 1e-6, reltol = 1e-6)
sols = sol[C_c_1(t, r)]
plot(sols[])