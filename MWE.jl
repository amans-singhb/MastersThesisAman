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


### PDE system ###

D_cat_val = 0.25e-3 # [m]
rad_cat = 0.5 * D_cat_val # [m]
L_val = 4.8e-3 # [m]

T_val = 450.0 # [K]
T_c_val = 500.0 # [K]
C_i_val = [5.0, 7.0]
drTc_val = 5

C_c_i_init = [0.1, 33.39396610065385]

## Parameters ##
@parameters t r T T_c drTc C_i[1:2] C_c_i[1:2]

## Differential ##
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_c_i(..)[1:2] T(..)[1:2] X(..) Y(..)

## Equations and Differential Equations ##

eqs = [Dt(C_c_i(t, r)[i]) ~ 1 * Dr(C_c_i(t, r)[i]) + Drr(C_c_i(t, r)[i]) for i in 1:2]

ICS_C_c_i = [C_c_i(0.0, r)[i] ~ C_c_i_init[i] for i in 1:2]

BCS2 = [Dr(C_c_i(t, 0)[i]) ~ 0 for i in 1:2]
BCS3 = [(C_c_i(t, rad_cat)[i] - C_i[i]) ~ Dr(C_c_i(t, rad_cat)[i]) for i in 1:2]
STEP1_BCS4 = [(C_c_i(t, rad_cat)[i] - C_i[i]) for i in 1:2]
STEP2_BCS4 = [Dr(C_c_i(t, rad_cat)[i]) for i in 1:2] # the sum won't work without a for loop
BCS4 = [(T_c - T) + sum(STEP1_BCS4) ~ drTc - sum(STEP2_BCS4)]

bcs = [ICS_C_c_i; BCS2; BCS3; BCS4]

using OrdinaryDiffEq, DomainSets, MethodOfLines
using ModelingToolkit: scalarize

params_vec = [C_i[i] => C_i_val[i] for i in 1:2]

# Domain
domains = [t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, rad_cat)]

# System
vars1 = reduce(vcat, [[C_c_i(t, r)[i] for i in 1:2], X(t, r)])
vars2 = reduce(vcat, [[T(t, r)[i] for i in 1:2], Y(t, r)])
vars = [vars1; vars2]
params_scal = [T => T_val, T_c => T_c_val, drTc => drTc_val]
params = [params_scal; params_vec]
@named WGS_pde = PDESystem(eqs, bcs, domains, [t, r], vars, params)

# Discretization
dr = round(rad_cat/10, sigdigits=4)
order = 2
discretization = MOLFiniteDifference([r => dr], t, order = order)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)
sol = solve(prob, Tsit5())