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

## Heat capacity ##

# Heat capacity for i [J/mol K]
function C_p_i_func(T)
    # p = [i, A, B, C, D]
    p = [1 6.60 1.20e-3 0 0;
        2 10.34 2.74e-3 0 -1.955e5;
        3 6.62 0.81e-3 0 0;
        4 8.22 0.15e-3 1.34e-6 0;
        5 6.50 1.00e-3 0 0]

    A, B, C, D = p[:, 2], p[:, 3], p[:, 4], p[:, 5]

    C_p_i = A + B * T + C * T^2 + D / T^2
    
    C_p_i = C_p_i * 4.184 # for conversion to [J/mol K]

    return C_p_i
end

# Heat capacity of mixture [J/mol K]
function C_p_func(C_i, T)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    
    C_p_i = C_p_i_func(T)
    
    C_p = sum(y .* C_p_i)

    return C_p
end

# Heat capacity of mixture [kcal/kmol K]
function C_p_cal_func(C_i, T)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    
    C_p_i = C_p_i_func(T) * (1/4.184)
    
    C_pc = sum(y .* C_p_i)

    return C_pc
end

# Heat capacity of catalyst [J/mol K]
function C_p_cat_func(T)
    C_p_cat = (0.435157281 * (7.39634088 + 0.00328222 * T - 300420.37 / T^2 + 88.2034012 / (T)^0.5) + 0.475928864 * (9.27078814 + 0.01612958 * T - 5.319e-05*T^2 + 1.0342e-07 * T^3) + 0.088913855 * ( 17.5417116 + 0.05947003 * T - 0.0001893*T^2 + 1.6385e-07 * T^3 + 2.6593e-09 * T^4)) # [cal/mol K]

    C_p_cat = C_p_cat * 4.184 # for conversion to [J/mol K]

    return C_p_cat
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

## Thermal conductivity ##

# Thermal conductivity of gas phase of species i [kcal/h m K]
function λ_i_func(T)
    # p = [i, A, B, C, D]
    p = [1 5.1489e-4 0.6863 57.13 501.92;
        2 3.1728 -0.3838 964 1.86e6;
        3 2.2811e-3 0.7452 12 0;
        4 5.3345e-6 1.3973 0 0;
        5 2.8497e-4 0.7722 16.323 373.72]

    A, B, C, D = p[:, 2], p[:, 3], p[:, 4], p[:, 5]

    λ_i =  @. (A * T^B) / (1 + C / T + D / T^2)

    return λ_i
end

# Binary interaction parameter A_ij 
function A_ij_func(i, j, μ_i, M_i, T, T_boil)
    C = 1.0 # [-]
    if i == j
        return 1.0
    else
        S_i = 1.5 * T_boil[i]
        S_j = 1.5 * T_boil[j]
        S_ij = C * (S_i * S_j)^0.5

        A_ij = (1/4) * (1 + ((μ_i[i]/μ_i[j]) * (M_i[j]/M_i[i])^0.75 * ((T + S_i) / (T + S_j)))^0.5 )^2 * ((T + S_ij) / (T + S_i))
        return A_ij
    end
end

# Thermal conductivity of mixture at atmospheric pressure [kcal/h m K] (also works above atmospheric pressure)
function λ_dash_func(C_i, T)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    
    M_i = [28.01, 44.01, 2.016, 18.016, 28.014] # [kg/kmol]
    T_boil = [81.65, 194.7, 20.35, 373, 77.36] # [K]
    μ_i = μ_i_func(T) # [kg/h m]
    λ_i = λ_i_func(T) # [kcal/h m K]

    λ_dash = 0

    for i in eachindex(y)
        i = i[1]
        numerator = y[i] * λ_i[i]
        denominator = 0
        for j in eachindex(y)
            j = j[1]
            denominator += y[j] * A_ij_func(i, j, μ_i, M_i, T, T_boil)
        end
        λ_dash += numerator / denominator
    end

    return λ_dash 
end

## Enthalpy ##

# Heat capacity integrated for i [J/mol] (enthalpy for i at T)
function C_p_i_integrated_func(T)
    # p = [i, A, B, C, D]
    p = [1 6.60 1.20e-3 0 0;
        2 10.34 2.74e-3 0 -1.955e5;
        3 6.62 0.81e-3 0 0;
        4 8.22 0.15e-3 1.34e-6 0;
        5 6.50 1.00e-3 0 0]

    A, B, C, D = p[:, 2], p[:, 3], p[:, 4], p[:, 5]

    C_p_i_integrated = A * T + B/2 * T^2 + C/3 * T^3 - D / T

    C_p_i_integrated = C_p_i_integrated * 4.184 # for conversion to [J/mol]

    return C_p_i_integrated
end

# Enthalpy of i [J/mol]
function H_i_func(T)
    H_form = [-110.53, -393.51, 0, -241.83, 0] * 1000 # for conversion to [J/mol]

    C_p_i_integrated_298 = C_p_i_integrated_func(298)
    C_p_i_integrated_T = C_p_i_integrated_func(T)

    H_i = H_form + (C_p_i_integrated_T - C_p_i_integrated_298)

    return H_i
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
    M_i = [28.01, 44.01, 2.016, 18.016, 28.014] # [kg/kmol]
    M = sum(y .* M_i)

    ρ = (P * M) / (R * T) # [kg/m^3]

    μ = μ_mix_func(y, T, M_i) # [kg/h m]
    
    ϵ_b = ϵ_b_func(D_rct, D_cat)
    G = G_func(F, D_rct, ϵ_b, M) # [kg/m^2 h]

    # k_c_i_1 =  0.357 * (((ρ * D_i_m) / μ)).^(2/3) * (G / (ρ * ϵ_b)) * ((μ / (D_cat * G)))^0.359 # in paper
    k_c_i_2 = 0.357 * G * ((ρ * D_i_m) / μ).^(2/3) / (ρ * ϵ_b * ((D_cat * G / μ)^0.359)) # in Jacobian files code (Very small difference in results, 6e-14 difference in k_c_i[1])

    return k_c_i_2
end

# Heat transfer coefficient [J/h m^2 K]
function h_f_func(T, C_i, D_cat, D_rct, F)
    y = [C_i[1]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[2]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[3]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[4]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);
        C_i[5]/(C_i[1] + C_i[2] + C_i[3] + C_i[4] + C_i[5]);]
    M_i = [28.01, 44.01, 2.016, 18.016, 28.014] # [kg/kmol]
    M = sum(y .* M_i)

    Cpc = C_p_cal_func(C_i, T)

    μ = μ_mix_func(y, T, M_i) # [kg/h m]

    λ = λ_dash_func(C_i, T) # [kcal/h m K]

    ϵ_b = ϵ_b_func(D_rct, D_cat)
    G = G_func(F, D_rct, ϵ_b, M) # [kg/m^2 h]

    # h_f_1 = @. 1.37 * (0.357 / ϵ_b) * ((Cpc * G) / M) * (μ / ( D_cat * G))^0.359 * ((λ * M)/ (Cpc * μ))^(2/3) # in the paper [kcal/h m^2 K]
    h_f_2 = @. 1.37 * 0.357 * (Cpc * G / M) / (ϵ_b * (D_cat * G / μ)^0.359 * ((Cpc * μ) / (λ * M))^(2/3) ) # in Jacobian files code [kcal/h m^2 K]

    h_f_2 = h_f_2 * 4184 # for conversion to [J/h m^2 K]

    return h_f_2
end

using DelimitedFiles

# util function to write data to csv
function write_to_csv(filename::String, data, folder_path::String, delimiter::String=",")
    if !isdir(folder_path)
        mkdir(folder_path)
    end
    file_path = joinpath(folder_path, filename)
    writedlm(file_path, data, delimiter)
end