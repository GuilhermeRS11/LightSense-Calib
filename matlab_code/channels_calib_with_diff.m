clc;        
clear all;  
close all;  

% Initialization of calibration values based on standard OSRAM data
F1 = 0.011057; F2 = 0.019664; F3 = 0.028302; F4 = 0.032492;
F5 = 0.034778; F6 = 0.034016; F7 = 0.039658; F8 = 0.043168;
Clear = 0.127734; NIR = 0.051806;

% Global variable storage for use in functions
global sensor_values matrix_GSCM;
sensor_values = [F1, F2, F3, F4, F5, F6, F7, F8, NIR];

% Configuration of Wavelengths and Correction Matrix

WaveLenght = 380:1:1000; % Total available wavelengths

% Loading and configuration of the spectral correction matrix
matrix_GSCM = xlsread("General_Spec_Corr_Matrix.xlsx");
matrix_GSCM = matrix_GSCM(2:end, 2:end);

% Application of gain corrections to adjust sensor values
gain_correction_256 = 0.987308373187068; % Gain for 256
gain_correction_512 = 0.959349243600411; % Gain for 512
gain_correction = [gain_correction_512(ones(1,6)), gain_correction_256(ones(1,4))];

% Load the average coefficients defined in the channel calibration
load('calibr_coeficients.mat');

%% Adjustment of Basic Counters Using the Reconstructed Spectrum

% Define the trials that will be analyzed
experiment_name = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", ...
                 "GW_100", "RGB_100", "RGBW_100"];

% Configuration for the number of objectives (types of error) and diffusers
n_opt_func = 4; % MSE, WMSE, MAE, RMSE
n_diffs = 2;            % Diffuser 1 and Diffuser 2
n_sersor_channels = 10; % Number of sensor channels

% Preparation to store the mean coefficients obtained by calibration
mean_coeficients = zeros(n_opt_func, n_diffs, n_sersor_channels);

% Optimization options settings for the lsqnonlin algorithm
options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt',... 
    'OptimalityTolerance', 1e-15, 'FunctionTolerance', 1e-50, 'MaxIterations', 1000, ...
    'StepTolerance', 1e-15, 'MaxFunctionEvaluations', 2.000000e+04);

% Iteration over all trials
for i = 1:length(experiment_name)
    for j = 1:n_diffs
        % Load sensor data with diffuser
        load(append("calibration_tests/sensor_measurements", experiment_name(i), sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_with_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
        
        % Load data without diffuser for comparison
        load(append("calibration_tests/sensor_measurements", experiment_name(i), "-100_amostras.mat"));
        sensor_values_without_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
        
        % Original data obtained from the combination of correction matrix and sensor values without diffuser
        original_data = matrix_GSCM * sensor_values_without_diffuser';
        
        % Function that reconstructs the spectrum from the adjusted parameters
        funcaoSPD = @(params) spectrum_restore(params, sensor_values_with_diffuser);
        
        % Iteration over each objective function
        for k = 1:n_opt_func
            % Selection of the objective function based on the type of error
            switch k
                case 1 % MSE
                    objective_function = @(params) mean((funcaoSPD(params) - original_data).^2);
                case 2 % WMSE
                    objective_function = @(params) mean((funcaoSPD(params) - original_data).^2 .* original_data');
                case 3 % MAE
                    objective_function = @(params) mean(abs(funcaoSPD(params) - original_data));
                case 4 % RMSE
                    objective_function = @(params) sqrt(mean((funcaoSPD(params) - original_data).^2));
            end
            
            % Initial guess and parameter adjustment
            initial_parameters = ones(1, n_sersor_channels);
            adjusted_parameters = lsqnonlin(objective_function, initial_parameters, [], [], options);
            
            % Storage of adjusted coefficients
            mean_coeficients(k, j, :) = adjusted_parameters;
        end
    end
end

%% Calculation of R² for each trial using each of the methods

% Definition of trials and test configurations
experiment_name = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", "GW_100", "RGB_100", "RGBW_100"];
experiment_labels = ["R", "G", "B", "W", "RB", "RG", "GW", "RGB", "RGBW"];
num_experiments = numel(experiment_name);
n_opt_func = 5;       % MSE, WMSE, MAE, RMSE, without calibration
n_difusores = 2;      % Diffuser 1 and Diffuser 2
R2_results = zeros(num_experiments, n_opt_func, n_difusores); % Matrix to store R² results

% Loop to process each trial
for i = 1:num_experiments
    experiment = experiment_name(i);
    
    % Load sensor data without diffuser for comparison
    load(append("calibration_tests/sensor_measurements", experiment, "-100_amostras.mat"));
    sensor_values_without_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
    
    % Obtain the reconstructed spectrum without using a diffuser for reference
    reconstructed_semDifo = spectrum_restore(ones(1,10), sensor_values_without_diffuser);
    y_mean = mean(reconstructed_semDifo);
    SS_tot = sum((reconstructed_semDifo - y_mean).^2);

    % Processing for each diffuser configuration
    for j = 1:n_difusores
        load(append("calibration_tests/sensor_measurements", experiment, sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_with_diffuser = sum(buffer_sensor_values,1)./99 .* gain_correction .* calibr_coeficientes;

        % Calculation of R² for the case without calibration adjustment
        reconstructed = spectrum_restore(ones(1,10), sensor_values_with_diffuser);
        SS_res = sum((reconstructed_semDifo - reconstructed).^2);
        R2_results(i, 1, j) = 1 - SS_res / SS_tot;

        % Calculation of R² for each calibration method
        for k = 2:n_opt_func
            mean_coef = squeeze(mean_coeficients(k-1, j, :))'; % Adjustment to align with matrices
            reconstructed = spectrum_restore(mean_coef, sensor_values_with_diffuser);
            SS_res = sum((reconstructed_semDifo - reconstructed).^2);
            R2_results(i, k, j) = 1 - SS_res / SS_tot;
        end
    end
end

%% Function that reconstructs the spectrum through the Golden Device matrix
function SPD_reconstructed = spectrum_restore(calibration_coeficients, sensor_values)
    % Import necessary global variables for the function
    global matrix_GSCM;
    
    % Definition of limits for calibration coefficients
    lower_limits = ones(1, 10) * -1;
    upper_limits = ones(1, 10) * 10;
    
    % Check if the calibration coefficients are within the established limits
    if any(calibration_coeficients < lower_limits) || any(calibration_coeficients > upper_limits)
        disp("Coefficient provided is out of allowed limits");
        % Return an error vector if the coefficients are out of limits
        SPD_reconstructed_norm = ones(size(WaveLenght)) * 1e10;
        return;
    end

    % If the coefficients are within the limits, perform the calculation
    % Multiply sensor values by calibration coefficients
    sensor_values_calibr = sensor_values .* calibration_coeficients;

    % Multiply the spectral correction matrix by the calibrated sensor values
    SPD_reconstructed = matrix_GSCM * sensor_values_calibr';

    % Normalize the reconstructed spectrum by its maximum value
    SPD_reconstructed_norm = SPD_reconstructed / max(SPD_reconstructed);
end
