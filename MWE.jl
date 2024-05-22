#### MWE for WGS reactor model ####

### Define functions ###

using Symbolics

# Binary gas diffusivities for component pairs {T [K], P [atm], D_ij [cm^2/s]}
function D_ij_func_a(T, P, A, B, C, D, E, F)
    (A * abs(T)^B / P) * abs((log(C / abs(T))))^(-2 * D) * exp(-E / T - F / abs(T)^2)
end

function D_ij_func_b(P, B)
    B / P
end

function D_ij_func_c(T, P, A, B)
    (A * T + B) / P
end

# Matrix of all binary gas diffusivities for component pairs
function D_ij_matrix_func(T, P, D_ij_matrix) 
    #p = [eq, i, j, A,      B,      C,  D,  E,   F] (C is set to 1 where it has no value, to avoid log(0) error, D accounts for lack of C in eq = "a" )
    p = ["a" 3 1 15.39e-3 1.548 0.316e8 1 -2.80 1067;
        "a" 3 2 3.14e-5 1.75 1 0 11.7 0;
        "b" 3 4 0 1.020 1 0 0 0;
        "c" 3 5 6.007e-3 -0.99311 1 0 0 0;
        "a" 1 2 3.15e-5 1.57 1 0 113.6 0;
        "a" 1 4 0.187e-5 2.072 1 0 0 0;
        "c" 1 5 0 0.322 1 0 0 0;
        "a" 2 4 9.24e-5 1.5 1 0 307.9 0;
        "a" 2 5 3.15e-5 1.57 1 0 113.6 0;
        "a" 4 5 0.187e-5 2.072 1 0 0 0;]
    
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

# Effective diffusivity of i in j
function D_eff_ij_func(D_ij, θ, τ)
    D_ij .* (θ / τ)
end

# Effective diffusivity of i in the mixture
function D_i_m_func(y, θ, τ, D_i_m_vec, T, P, D_ij_matrix)
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

    return D_i_m_vec
end

# Viscosity of gas phase of species i {T [K], μ [N s/m^2]}
function μ_i_func(T, i)
    # p = [i, A, B, C, D]
    p = [1 1.1127e-6 0.5338 94.7 0;
        2 2.148e-6 0.46 290 0;
        3 1.797e-7 0.685 -0.59 140;
        4 1.7096e-8 1.1146 0 0;
        5 6.5592e-7 0.6081 54.714 0]

    row = p[i, :]
    A, B, C, D = row[2], row[3], row[4], row[5]

    return (A * abs(T)^B) / (1 + C / T + D / abs(T)^2)
end

# Viscosity of gas phase mixture
function μ_mix_func(y, μ_i, M_i)
    μ_mix = 0

    for i in eachindex(y)
        i = i[1]
        numerator = y[i] * μ_i[i]
        denominator = 0
        for j in eachindex(y)
            j = j[1]
            denominator += y[j] * sqrt(M_i[j] / M_i[i])
        end
        μ_mix += numerator / denominator
    end

    return μ_mix
end

# Reynolds number
function Re_func(ρ, u, L, μ)
    (ρ * u * L) / μ
end

# Friction factor correlation for gas flow through packed tubular reactor
function F_fr_func(G, D_cat, ρ, ϵ_b, Re)
    (- G^2 / (ρ * D_cat)) * ((1 - ϵ_b) / ϵ_b^3) * (1.75 + 4.2 * ((1 - ϵ_b) / abs(Re)^(1/6)))   
end

# Molar density of the gas phase
function ρ_func(P, T, R)
    P / (R * T)
end

# Molar fraction of i in the gas phase
function y_func(C_i, i)
    C_i[i] / sum(C_i)
end

# Linear gas velocity
function u_func(α, T, P)
    α * (T / P)
end


# Molar flux
function G_func(F_0, D_rct, ϵ_b)
    (4 * F_0) / (pi * D_rct^2 * ϵ_b)
end

# Alpha (momentum balance term)
function α_func(G, R)
    G * R
end

# Bed porosity
function ϵ_b_func(D_rct, D_cat)
    0.38 + 0.073 * (1 - (D_rct / D_cat - 2)^2 / (D_rct / D_cat)^2)
end


### Register symbolic functions ###
using ModelingToolkit

@register_symbolic D_ij_func_a(T, P, A, B, C, D, E, F)
@register_symbolic D_ij_func_b(P, B)
@register_symbolic D_ij_func_c(T, P, A, B)
@register_symbolic D_ij_matrix_func(T, P, D_ij_matrix)
@register_symbolic D_eff_ij_func(D_ij, θ, τ)
@register_symbolic D_i_m_func(y, θ, τ, D_i_m_vec, T, P, D_ij_matrix)

@register_symbolic μ_i_func(T, i)
@register_symbolic μ_mix_func(y, μ_i, M_i)

@register_symbolic Re_func(ρ, u, L, μ)
@register_symbolic F_fr_func(G, D_cat, ρ, ϵ_b, Re)
@register_symbolic ρ_func(P, T, R)
@register_symbolic y_func(C_i, i)
@register_symbolic u_func(α, T, P)


### PDE system ###

# Parameters for inlet and init
F_0 = 10 # [mol/h]

R_val = 8.2057e-5 # [m3 atm/mol K]

y_0 = [0.2749, 0.1198, 0.2686, 0.3165, 0.0202]
y_init = [0.0, 0.0, 0.0, 0.0, 1.0]


# Inlet conditions
T_in = 473 # [K] 200 C
P_in = 1.3 # [atm]
V_flow_0 = (F_0 * R_val * T_in) / P_in # [m3/h]
C_i_in = y_0 * (F_0 / V_flow_0) # [mol/m3]

# Initial conditions
T_init = 473 # [K] 200 C
P_init = 1.3 # [atm]
V_flow_init = (F_0 * R_val * T_init) / P_init # [m3/h]
C_i_init = y_init * (F_0 / V_flow_init) # [mol/m3]

# Other parameters
M_i_val = [28.01, 44.01, 2.016, 18.016, 28.014]
θ_val = 0.55 #[-]
τ_val = 5 #[-]

D_cat_val = 0.25e-3 # [m]
rad_cat = 0.5 * D_cat_val # [m]
D_rct = 12.7e-3 # [m]
L_val = 4.8e-3 # [m]

ϵ_b_val = ϵ_b_func(D_rct, D_cat_val)
G_val = G_func(F_0, D_rct, ϵ_b_val)
α_val = α_func(G_val, R_val)

## Parameters ##
@parameters t z r θ τ α R M_i[1:5] G D_cat ϵ_b L

## Differential ##
Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
Drr = Differential(r)^2

## Variables ##
@variables C_i(..)[1:5], T(..)[1:5], P(..)[1:5], C_c_i(..)[1:5], T_c(..)[1:5], y_c(t, z, r)[1:5], y(t, z)[1:5], D_i_m(t, z, r)[1:5], ρ(t, z), Re(t, z), μ_i(t, z)[1:5], μ(t, z), u(t, z)

## Equations and Differential Equations ##
import ModelingToolkit: scalarize

equations_y = [y[i] ~ y_func(C_i(t, z), i) for i in 1:5]
equations_y_c = [y_c[i] ~ y_func(C_c_i(t, z, r), i) for i in 1:5]
equations_μ = [μ_i[i] ~ μ_i_func(T(t, z)[i], i) for i in 1:5]

equations3 = [ρ ~ ρ_func(P(t, z)[1], T(t, z)[1], R),
    μ ~ μ_mix_func(y, μ_i, M_i),
    u ~ u_func(α, T(t, z)[1], P(t, z)[1]),
    Re ~ Re_func(ρ, u, L, μ)]

matrix_D_ij = zeros(Num, 5, 5)
vector_D_i_m = zeros(Num, 5)
equations_D_i_m = [scalarize(D_i_m[1:5] .~ D_i_m_func(y_c, θ, τ, vector_D_i_m, T_c(t, z, r)[1], P(t, z)[1], matrix_D_ij))...]

equations = [equations_y_c; equations_y; equations_μ; equations3; equations_D_i_m]

DE1 = [Dt(C_i(t, z)[i]) ~ Dz(C_i(t, z)[i]) - (Dz(T(t, z)[i]) - Dz(P(t, z)[i])) for i in 1:5]
DE2 = [Dz(P(t, z)[i]) ~ - F_fr_func(G, D_cat, ρ, ϵ_b, Re) for i in 1:5]
DE3 = [Dt(T(t, z)[i]) ~ Dz(T(t, z)[i]) for i in 1:5]
DE4 = [Dt(C_c_i(t, z, r)[i]) ~ (Dr(D_i_m[i])) * Dr(C_c_i(t, z, r)[i]) + Drr(C_c_i(t, z, r)[i]) for i in 1:5]
STEP_DE5 = [Dt(T_c(t, z, r)[i]) for i in 1:5]
DE5 = [Dt(T_c(t, z, r)[i]) + sum(STEP_DE5) ~ (Dr(T_c(t, z, r)[i]) + Drr(T_c(t, z, r)[i])) + (Dr(C_c_i(t, z, r)[i]) * Dr(T_c(t, z, r)[i])) for i in 1:5]

diffequations = [DE1; DE2; DE3; DE4; DE5]

eqs = [equations; diffequations]

## Initial and Boundary Conditions ##
ICS_T = [T(0.0, z)[i] ~ T_init for i in 1:5]
ICS_P = [P(0.0, z)[i] ~ P_init for i in 1:5]
ICS_C_i = [C_i(0.0, z)[i] ~ C_i_init[i] for i in 1:5]

ics = [ICS_T; ICS_P; ICS_C_i]

BCS_T = [T(t, 0.0)[i] ~ T_in for i in 1:5]
BCS_P = [P(t, 0.0)[i] ~ P_in for i in 1:5]
BCS_Tc = [Dz(T_c(t, z, 0)[i]) ~ 0 for i in 1:5]
BCS1 = [C_i(t, 0.0)[i] ~ C_i_in[i] for i in 1:5]
BCS2 = [Dz(C_c_i(t, z, 0)[i]) ~ 0 for i in 1:5]
BCS3 = [(C_c_i(t, z, rad_cat)[i] - C_i(t, z)[i]) ~ Dr(C_c_i(t, z, rad_cat)[i]) for i in 1:5]
STEP_BCS4 = [Dr(C_c_i(t, z, rad_cat)[i]) for i in 1:5] # the sum won't work without a for loop
BCS4 = [(T_c(t, z, rad_cat)[i] - T(t, z)[i]) + sum((C_c_i(t, z, rad_cat) - C_i(t, z))) ~ Dr(T_c(t, z, rad_cat)[i]) - sum(STEP_BCS4) for i in 1:5]

bcs = [ics; BCS_T; BCS_P; BCS_Tc; BCS1; BCS2; BCS3; BCS4]

using OrdinaryDiffEq, DomainSets, MethodOfLines

# Domain
domains = [t ∈ Interval(0.0, 1.0),
    z ∈ Interval(0.0, L_val),
    r ∈ Interval(0.0, rad_cat)]

# System
vars = [C_i(t, z)[1:5], T(t, z)[1:5], P(t, z)[1:5], C_c_i(t, z, r)[1:5], T_c(t, z, r)[1:5], y_c[1:5], y[1:5], D_i_m[1:5], ρ, Re, μ_i[1:5], μ, u]
params = [θ => θ_val, τ => τ_val, α => α_val, R => R_val, M_i[1:5] => M_i_val[1:5], G => G_val, D_cat => D_cat_val, ϵ_b => ϵ_b_val, L => L_val]

@named WGS_pde = PDESystem(eqs, bcs, domains, [t, z, r], vars, params)

# Discretization
dz = round(L_val/10, sigdigits=4)
dr = round(rad_cat/10, sigdigits=4)
order = 2
discretization = MOLFiniteDifference([z => dz, r => dr], t)

# Converting PDE to ODE with MOL
prob = discretize(WGS_pde, discretization)