#### Functions for catalyst particle balance for WGS reactor model ####

using Pkg
Pkg.activate("WGS")

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

### Define functions for physical properties ###

using Symbolics

## Diffusivity ##

# Binary gas diffusivities for component pairs [cm^2/s]
function D_ij_func_a(T, P, A, B, C, D, E, F)
    (A * T^B / P) * (log(C / T))^(-2 * D) * exp(-E / T - F / T^2)
end

function D_ij_func_b(P, B)
    B / P
end

function D_ij_func_c(T, P, A, B)
    (A * T + B) / P
end

# Matrix of all binary gas diffusivities for component pairs [cm^2/s]
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

# Effective diffusivity of i in the mixture [m^2/h]
function D_i_m_func(C_i, θ, τ, T, P)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    
    D_ij_matrix = zeros(Num, 5, 5)
    D_i_m_vec = zeros(Num, 5)

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

    D_i_m_vec = D_i_m_vec * (3600/10000) # for conversion to [m^2/h]

    return D_i_m_vec
end

## Viscosity ##

# Viscosity of gas phase of species i [N s/m^2]
function μ_i_func(T)
    # p = [i, A, B, C, D]
    p = [1 1.1127e-6 0.5338 94.7 0;
        2 2.148e-6 0.46 290 0;
        3 1.797e-7 0.685 -0.59 140;
        4 1.7096e-8 1.1146 0 0;
        5 6.5592e-7 0.6081 54.714 0]

    A, B, C, D = p[:, 2], p[:, 3], p[:, 4], p[:, 5]

    μ_i = @. (A * T^B) / (1 + C / T + D / T^2)
    
    return μ_i
end

# Viscosity of gas phase mixture [kg/h m]
function μ_mix_func(y, T, M_i)
    μ_i = μ_i_func(T)

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

    μ_mix = μ_mix * 3600 # for conversion to [kg/h m]

    return μ_mix
end

## Other functions ##

# Reaction rate of i [mol/m^3 h]
function r_i_func(C_i, d_cat, θ, P, T)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]

    
    K_eq = exp((4577.8 / T) - 4.33)
    R_joule = 8.314
    F_press = P^(0.5-P/500)
    r_co_min = (-1) * d_cat * (1 - θ) * F_press * (2.96e5) * exp(-47400 / (R_joule * T)) * P^2 * (y[1] * y[4] - ((y[2] * y[3]) / K_eq))
    r_co_min = r_co_min * 1000 # for conversion to [mol/m^3 h]
    
    #       CO,         CO2,        H2,     H2O,    N2
    r_i = [r_co_min, -r_co_min, -r_co_min, r_co_min, 0]

    return r_i
end

# Bed porosity
function ϵ_b_func(D_rct, D_cat)
    0.38 + 0.073 * (1 - ((D_rct / D_cat - 2)^2) / ((D_rct / D_cat)^2))
end

# Molar flux [kg/m^2 h]
function G_func(F_0, D_rct, ϵ_b, M)
    (0.001 * 4 * F_0 * M) / (pi * D_rct^2 * ϵ_b) # 0.001 for conversion to [kg/m^2 h]
end

# Mass transfer coefficient [m/h]
function k_c_i_func(T, P, R, C_i, D_i_m, D_cat, D_rct, F)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    M_i = [28.01, 44.01, 2.016, 18.016, 28.014]
    M = sum(y .* M_i)

    ρ = (P * M) / (R * T) # [kg/m^3]

    μ = μ_mix_func(y, T, M_i) # [kg/h m]
    
    ϵ_b = ϵ_b_func(D_rct, D_cat)
    G = G_func(F, D_rct, ϵ_b, M) # [kg/m^2 h]

    # k_c_i_1 =  0.357 * (((ρ * D_i_m) / μ)).^(2/3) * (G / (ρ * ϵ_b)) * ((μ / (D_cat * G)))^0.359 # in paper
    k_c_i_2 = 0.357 * G * ((ρ * D_i_m) / μ).^(2/3) / (ρ * ϵ_b * ((D_cat * G / μ)^0.359)) # in Jacobian files code (Very small difference in results, 6e-14 difference in k_c_i[1])

    return k_c_i_2
end

## Util functions ##

using DelimitedFiles

# util function to write data to csv
function write_to_csv(filename::String, data, folder_path::String, delimiter::String=",")
    if !isdir(folder_path)
        mkdir(folder_path)
    end
    file_path = joinpath(folder_path, filename)
    writedlm(file_path, data, delimiter)
end

# Function for extracting data from solutions
function extractData(sol; t_length = length(sol.t), r_length = 21, modelPDESize = 5, idx = [10, 4, 3, 11, 6])
    data = [zeros(Float64, t_length, modelPDESize) for _ in 1:r_length]
    t_range = sol.t
    r_range = range(0.0, 0.125e-3, length=r_length)
    # i for radial position, j for time
    for i in 1:r_length
        for j in 1:t_length
            data[i][j,:] = sol(t_range[j], r_range[i])[idx]
        end
    end
    return data
end

# function for numerical differentiation of solutions (forward, backward and central difference)
function finite_diff_sol(sol)
    sol_cci = [sol[C_c_1(t, r)], sol[C_c_2(t, r)], sol[C_c_3(t, r)], sol[C_c_4(t, r)], sol[C_c_5(t, r)]]
    sol_time = sol.t

    num_time = length(sol_time)
    num_comp = length(sol_cci)
    size_sol = size(sol_cci[1])

    Δcci_t = [zeros(size_sol) for _ in 1:num_comp]
    err = [zeros(num_time) for _ in 1:num_comp]

    for i in 1:num_comp
        Δcci_t[i][1, :] = (sol_cci[i][2, :] - sol_cci[i][1, :]) / (sol_time[2] - sol_time[1])
        err[i][1] = sol_time[2] - sol_time[1]

        for j in 2:(num_time-1)
            Δcci_t[i][j, :] = (sol_cci[i][j+1, :] - sol_cci[i][j-1, :]) / (sol_time[j+1] - sol_time[j-1])
            err[i][j] = (sol_time[j+1] - sol_time[j-1])^2
        end

        Δcci_t[i][end, :] = (sol_cci[i][end, :] - sol_cci[i][end-1, :]) / (sol_time[end] - sol_time[end-1])
        err[i][end] = sol_time[end] - sol_time[end-1]
    end

    return Δcci_t, err
end