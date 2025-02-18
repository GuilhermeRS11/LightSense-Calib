clc;        
clear all;  
close all;  

% Initialization of calibration values based on standard OSRAM data
F1 = 0.011057; F2 = 0.019664; F3 = 0.028302; F4 = 0.032492;
F5 = 0.034778; F6 = 0.034016; F7 = 0.039658; F8 = 0.043168;
Clear = 0.127734; NIR = 0.051806;

% Global storage of variables for use in functions
global sensor_values WaveLenght matrix_GSCM;
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

% Defining the tests and initial configurations
experiment_name = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", "GW_100", "RGB_100", "RGBW_100"];
experiment_labels = ["R", "G", "B", "W", "RB", "RG", "GW", "RGB", "RGBW"];
num_experiments = numel(experiment_name);
n_opt_func = 2;  % lsqnonlin (1), particleswarm (2)
n_obj_func = 4;       % MSE, WMSE, MAE, RMSE
n_sersor_channels = 10;

%% Method 1 (Only 1 trial) - Adjust basic_counts through the reconstructed spectrum
% The commented code fragments can be used to analyze
% other objective functions and other optimization functions

experiment = 'RB_100';
global sphere_SPD;

% Load the sensor and integrating sphere test data
load(append("calibration_tests/sensor_measurements/", experiment,"-100_amostras.mat"));
sensor_values_no_diff = sum(buffer_sensor_values,1)./99;

% Apply gain calibration
sensor_values_no_diff = sensor_values_no_diff .* gain_correction;

% Adapt the sphere values from 800 to 1100 to be zero, as we don't have information for them
sphere_SPD = zeros(1, length(WaveLenght)); 
load(append("calibration_tests/sphere_measurements/spectrum/", experiment,"_S.mat"));
sphere_SPD(1:421) = data.Wrad_relative;
corr_coeff_no_diff = zeros(length(sensor_values_no_diff));

% Function that generates the reconstructed SPD from the parameters
SPD_function = @(params) spectrum_restore(params, sensor_values_no_diff);

% Objective functions to minimize the difference between the reconstructed SPD and the original SPD
%MSE (Mean Squared Error) unweighted
objective_function = @(params) mean((SPD_function(params) - sphere_SPD').^2); 

% WEQM (Weighted Mean Squared Error)
%objective_function = @(params) mean((SPD_function(params) - sphere_SPD').^2 .* sphere_SPD');

% MAE (Mean Absolute Error) 
%objective_function = @(params) mean(abs(SPD_function(params) - sphere_SPD'));

% RMSE (Root Mean Squared Error)
% objective_function = @(params) sqrt(mean((SPD_function(params) - sphere_SPD').^2));

% Initial guess for the 10 adjustable parameters
initial_parameters = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% Perform parameter adjustment using lsqnonlin
options = optimoptions('lsqnonlin','Display','iter'); % To display information about the optimization process
options.Algorithm = 'levenberg-marquardt';
options.OptimalityTolerance = 1e-15;
options.FunctionTolerance = 1e-50;
options.MaxIterations = 1000; % Set the maximum number of iterations
options.StepTolerance = 1e-15; % Set the step tolerance
options.MaxFunctionEvaluations = 1e4;
adjusted_parameters = lsqnonlin(objective_function, initial_parameters, [], [], options);

% lb = 0.01 * ones(1, 10);  % Lower bound for all coefficients
% ub = 10 * ones(1, 10);   % Upper bound for all coefficients
% options = optimoptions('particleswarm', 'SwarmSize', 30, 'HybridFcn', @fmincon, 'Display', 'iter');
% [adjusted_parameters, fval] = particleswarm(objective_function, 10, lb, ub, options);

SPD_reconstructed_original = SPD_function(initial_parameters);

% Use the adjusted parameters to obtain the optimized reconstructed SPD
SPD_reconstructed_adjusted = SPD_function(adjusted_parameters);
corr_coeff_no_diff = adjusted_parameters;      

MSE = mean((SPD_reconstructed_original - sphere_SPD').^2); % Mean Squared Error
WMSE = mean((SPD_reconstructed_original - sphere_SPD').^2.*sphere_SPD'); % Weighted Mean Squared Error
disp(append("MSE at the beginning = ",num2str(MSE)));
disp(append("WMSE at the beginning = ",num2str(WMSE)));

MSE = mean((SPD_reconstructed_adjusted - sphere_SPD').^2); % Mean Squared Error
WMSE = mean((SPD_reconstructed_adjusted - sphere_SPD').^2.*sphere_SPD'); % Weighted Mean Squared Error
disp(append("MSE at the end = ",num2str(MSE)));
disp(append("WMSE at the end = ",num2str(WMSE)));

% Plot results for evaluation
cmap = wavelengthToRGB(WaveLenght);

f = figure;
h1 = plot(WaveLenght, sphere_SPD,'--k','linew',2);
%title("Channel Calibration");
hold on;
h2 = plot(WaveLenght, SPD_reconstructed_original, 'r', 'LineWidth', 1);
h3 = plot(WaveLenght, SPD_reconstructed_adjusted, 'b', 'LineWidth', 1);

minValue = min(min(SPD_reconstructed_original), min(SPD_reconstructed_adjusted));
maxValue = 1;

% Add colored rectangle to represent the spectrum colors
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*(0.15), minValue - maxValue*(0.15), minValue - maxValue*(0.05), minValue - maxValue*(0.05)]; % Height of the rectangle
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

xlabel("Wavelength [nm]");
ylabel("Relative Magnitude");
grid on; 

leg = legend([h1, h2, h3], ["Sphere SPD", "Reconstruction", "Adjusted Reconstruction"], 'Location', 'northeast');
ylim([minValue - maxValue*(0.15), maxValue*(1.1)]);
xlim([WaveLenght(1) WaveLenght(end)]);

exportgraphics(f,append('saved_images/CalibrationIdealCoefficients_',experiment,'.png'),'ContentType','vector');

%% Method 1 (All trials) - Adjust basic_counts through the reconstructed spectrum

mean_coeficients = zeros(n_opt_func, n_obj_func, n_sersor_channels);
std_coeficients = zeros(n_opt_func, n_obj_func, n_sersor_channels);
corr_coeficients = zeros(n_opt_func, n_obj_func, num_experiments, n_sersor_channels);
tic;

% Process all data for the number of optimizer functions
for k = 1:n_opt_func
    % Process all data for the number of objective functions
    for j = 1:n_obj_func
        % Perform calibration for all requested experiments
        for i = 1:num_experiments  

            % Load the sensor and integrating sphere test data
            load(append("calibration_tests/sensor_measurements/", experiment_name(i),"-100_amostras.mat"));
            sensor_values_no_diff = sum(buffer_sensor_values,1)./99;

            % Apply gain calibration
            sensor_values_no_diff = sensor_values_no_diff .* gain_correction;

            % Adapt the sphere values from 800 to 1100 to be zero, as we don't have information for them
            sphere_SPD = zeros(1, length(WaveLenght)); 
            load(append("calibration_tests/sphere_measurements/spectrum/", experiment_name(i),"_S.mat"));
            sphere_SPD(1:421) = data.Wrad_relative;

            % Function that generates the reconstructed SPD from the parameters
            SPD_function = @(params) spectrum_restore(params, sensor_values_no_diff);
                     
            % Objective function to minimize the difference between the reconstructed SPD and the original SPD
            if j == 1
                % MSE (Mean Squared Error) unweighted
                objective_function = @(params) mean((SPD_function(params) - sphere_SPD').^2);
            
            elseif j == 2
                % WMSE (Weighted Mean Squared Error) 
                objective_function = @(params) mean(((SPD_function(params) - sphere_SPD').^2 ).* sphere_SPD');

            elseif j == 3
                % MAE (Mean Absolute Error) 
                objective_function = @(params) mean(abs(SPD_function(params) - sphere_SPD'));
                
            else    
                % RMSE (Root Mean Squared Error)
                objective_function = @(params) sqrt(mean((SPD_function(params) - sphere_SPD').^2));
            end
            
            % Initial guess for the 10 adjustable parameters
            initial_parameters = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

            if k == 1
                % Perform parameter optimization with lsqnonlin
                options = optimoptions('lsqnonlin','Display','iter'); % To display information about the optimization process
                options.Algorithm = 'levenberg-marquardt';
                options.OptimalityTolerance = 1e-15;
                options.FunctionTolerance = 1e-50;
                options.MaxIterations = 1000; % Set the maximum number of iterations
                options.StepTolerance = 1e-15; % Set the step tolerance
                options.MaxFunctionEvaluations = 1e4;
                adjusted_parameters = lsqnonlin(objective_function, initial_parameters, [], [], options);

            else
                % Perform parameter optimization with particleswarm
                lb = 0.01 * ones(1, 10);  % Lower bound for all coefficients
                ub = 10 * ones(1, 10);   % Upper bound for all coefficients
                options = optimoptions('particleswarm', 'SwarmSize', 30, 'HybridFcn', @fmincon, 'Display', 'iter');
                [adjusted_parameters, fval] = particleswarm(objective_function, 10, lb, ub, options);
            end
            
            corr_coeficients(k, j, i, :) = adjusted_parameters;  
        end
        % Calculate the mean along the dimension of the trials for each channel
        mean_coeficients(k, j, :) = mean(squeeze(corr_coeficients(k, j, :, :)), 1);
        std_coeficients(k, j, :) = std(squeeze(corr_coeficients(k, j, :, :)), 1);
    end
end
toc;

%% Calculation of R² for each experiment using each of the methods

% Number of experiments
num_experiments = numel(experiment_name);

% Correct initialization of variables
R2_MSE_values = zeros(n_opt_func, num_experiments);
R2_WMSE_values = zeros(n_opt_func, num_experiments);
R2_MAE_values = zeros(n_opt_func, num_experiments);
R2_RMSE_values = zeros(n_opt_func, num_experiments);
values_no_calib = zeros(1, num_experiments);
no_coeficients = ones(1, 10);

% Loop to calculate the metrics
for j = 1:n_opt_func
    for i = 1:num_experiments
        experiment = experiment_name(i);
        % Load sensor data
        load(append("calibration_tests/sensor_measurements/", experiment, "-100_amostras.mat"));
        sensor_values_no_diff = sum(buffer_sensor_values,1)./99;
        sensor_values_no_diff = sensor_values_no_diff .* gain_correction;
    
        % Load and adjust sphere data
        sphere_SPD = zeros(1, length(WaveLenght));
        load(append("calibration_tests/sphere_measurements/spectrum/", experiment, "_S.mat"));
        sphere_SPD(1:421) = data.Wrad_relative;
        
        % Metric to evaluate the performance of the methods through R²
        y_mean = mean(sphere_SPD);
        SS_tot = sum((sphere_SPD - y_mean).^2);
        
        % Reconstruction without calibration
        reconstructed_no_calib = spectrum_restore(no_coeficients, sensor_values_no_diff);
        SS_res_no_calib = sum((sphere_SPD' - reconstructed_no_calib).^2);
        values_no_calib(1,i) = 1 - SS_res_no_calib / SS_tot;
    
        % Reconstruct using coefficients and calculate metrics
        reconstructed = spectrum_restore(squeeze(mean_coeficients(j, 1, :))', sensor_values_no_diff);
        SS_res = sum((sphere_SPD' - reconstructed).^2);
        R2_MSE_values(j,i) = 1 - SS_res / SS_tot;
        
        reconstructed = spectrum_restore(squeeze(mean_coeficients(j, 2, :))', sensor_values_no_diff);
        SS_res = sum((sphere_SPD' - reconstructed).^2);
        R2_WMSE_values(j,i) = 1 - SS_res / SS_tot;
         
        reconstructed = spectrum_restore(squeeze(mean_coeficients(j, 3, :))', sensor_values_no_diff);
        SS_res = sum((sphere_SPD' - reconstructed).^2);
        R2_MAE_values(j,i) = 1 - SS_res / SS_tot;
         
        reconstructed = spectrum_restore(squeeze(mean_coeficients(j, 4, :))', sensor_values_no_diff);
        SS_res = sum((sphere_SPD' - reconstructed).^2);
        R2_RMSE_values(j,i) = 1 - SS_res / SS_tot;
    end
end

% Plotting the bars including spectrum without calibration
barWidth = 0.15;  % Adjustment to accommodate an extra bar
% Colors specified for each bar
barColors = [
    [0.3 0.3 0.3];  % Dark gray for 'Without calibration'
    [1 0.3 0.3];    % Soft red for MSE
    [0.5 1 0.5];    % Soft green for WMSE
    [1 0.5 1];      % Soft magenta for MAE
    [0.5 0.5 1];    % Soft blue for RMSE
];
% Reconstruct using coefficients from all methods for this experiment
positions = (1:num_experiments) + (i-1)*barWidth*4 - barWidth*(num_experiments/2-0.5); % Adjusted positions for the bars
R2_values_lsqnonlin = [R2_MSE_values(1,:) R2_WMSE_values(1,:) R2_MAE_values(1,:) R2_RMSE_values(1,:)];
R2_values_particleswarm = [R2_MSE_values(2,:) R2_WMSE_values(2,:) R2_MAE_values(2,:) R2_RMSE_values(2,:)];

titles_otim = ["lsqnonlin", "particleswarm"];

f = figure;
for j = 1:n_opt_func
    subplot(n_opt_func,1,j);
    hold on;
    bar(positions - 2 * barWidth, values_no_calib, barWidth, 'FaceColor',  barColors(1, :));  % Black for uncalibrated
    bar(positions - 1 * barWidth, R2_MSE_values(j,:), barWidth, 'FaceColor', barColors(2, :));
    bar(positions, R2_WMSE_values(j,:), barWidth, 'FaceColor', barColors(3, :));
    bar(positions + 1 * barWidth, R2_MAE_values(j,:), barWidth, 'FaceColor', barColors(4, :));
    bar(positions + 2 * barWidth, R2_RMSE_values(j,:), barWidth, 'FaceColor', barColors(5, :));
    
    % Graph settings
    set(gca, 'XTick', positions, 'XTickLabel', experiment_labels);
    xtickangle(45);
    legend({'Without calibration', 'MSE', 'WMSE', 'MAE', 'RMSE'}, 'Location', 'bestoutside');
    title([titles_otim(j)]);
    xlabel('Experiment');
    ylabel('R²');
    grid on;
    axis tight;
    min_value = min(min(R2_values_lsqnonlin), min(R2_values_particleswarm));
    ylim([min_value - abs(min_value * 0.5), 1.1]);
    hold off;
end

% Export the mean coefficients
calibr_coeficients = squeeze(mean_coeficients(1, 1, :))'; % Export the MSE that obtained the best result
save('calibr_coeficients.mat', 'calibr_coeficients');

set(gcf, 'Position', get(0, 'Screensize'));  % Maximizes the figure window
exportgraphics(f,'saved_images/All_R2.png','ContentType','vector');

%% Reconstructed spectrum using the average coefficient of each method

% Define the experiment to be analyzed. Within the standard of "Complete Experiments"
experiment = "RGBW_100";
Nome_grafico = "RGBW";
file = append("calibration_tests/sphere_measurements/spectrum/",experiment,"_S.txt")';

load(append("calibration_tests/sensor_measurements/",experiment,"-100_amostras.mat"));
sensor_values_no_diff = sum(buffer_sensor_values,1)./99;
sensor_values_no_diff = sensor_values_no_diff .* gain_correction;

% Adapt the sphere values from 800 to 1100 to be zero, as we don't have information for them
sphere_SPD = zeros(1, length(WaveLenght)); 
load(append("calibration_tests/sphere_measurements/spectrum/", experiment,"_S.mat"));
sphere_SPD(1:421) = data.Wrad_relative;


reconstruction_no_diff = matrix_GSCM * sensor_values_no_diff';
reconstruction_no_diff = reconstruction_no_diff ./ max(reconstruction_no_diff);

SPD_MSE = spectrum_restore(squeeze(mean_coeficients(1, 1, :))', sensor_values_no_diff);
SPD_WMSE = spectrum_restore(squeeze(mean_coeficients(1, 2, :))', sensor_values_no_diff);
SPD_MAE = spectrum_restore(squeeze(mean_coeficients(1, 3, :))', sensor_values_no_diff);
SPD_RMSE = spectrum_restore(squeeze(mean_coeficients(1, 4, :))', sensor_values_no_diff);

% To map the colors of the rectangle that shows the spectrum color
cmap = wavelengthToRGB(WaveLenght);

f = figure;
plot(WaveLenght, sphere_SPD,'--k','linew',2);
title(Nome_grafico);
hold on;

plot(WaveLenght, reconstruction_no_diff,'Color', barColors(1,:),'lineWidth', 1);
plot(WaveLenght, SPD_MSE,'Color', barColors(2,:),'lineWidth', 1);
plot(WaveLenght, SPD_WMSE,'Color', barColors(3,:),'lineWidth', 1);
plot(WaveLenght, SPD_MAE,'Color', barColors(4,:),'lineWidth', 1);
plot(WaveLenght, SPD_RMSE,'Color', barColors(5,:),'lineWidth', 1);


minValue = min([min(reconstruction_no_diff), min(SPD_MSE), min(SPD_WMSE), min(SPD_MAE), min(SPD_RMSE), min(sphere_SPD)]);
maxValue = max([max(reconstruction_no_diff), max(SPD_MSE), max(SPD_WMSE), max(SPD_MAE), max(SPD_RMSE), max(sphere_SPD)]);
% Add colored rectangle to represent the spectrum colors
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*(0.15), minValue - maxValue*(0.15), minValue - maxValue*(0.05), minValue - maxValue*(0.05)]; % Height of the rectangle (adjust as needed)
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

xlabel("Wavelength [nm]");
ylabel("Relative Magnitude");
grid on; 
legend(["Sphere SPD" "Uncalibrated" "MSE" "WMSE" "MAE" "RMSE"],'location','northeast');
ylim([minValue - maxValue*(0.15), maxValue*(1.1)]);
xlim([380, 1000]);

exportgraphics(f,append('saved_images/SPD_reconstructed_methods_',experiment,'.png'),'ContentType','vector');

%% Function that reconstructs the spectrum through the Golden Device matrix
function SPD_reconstruida_norm = spectrum_restore(calibration_coeficients, sensor_values)
    global WaveLenght matrix_GSCM;
    
    % Limits for calibration coefficients
    lower_limits = ones(1, 10) * 0.01;
    upper_limits = ones(1, 10) * 10;
    
    % Check if the coefficients are within the limits
    if any(calibration_coeficients < lower_limits) || any(calibration_coeficients > upper_limits)
        disp("Provided coefficient out of allowed limits");
        SPD_reconstruida_norm = ones(length(WaveLenght),1) * 1e10;  % Error output
        return;
    end

    % If the coefficients are within the limits, proceed with the calculation
    sensor_values_calibr = sensor_values .* calibration_coeficients;
    SPD_reconstruida = matrix_GSCM * sensor_values_calibr';
    SPD_reconstruida_norm = SPD_reconstruida / max(SPD_reconstruida);  % Normalize the result
      
end

% Function to map wavelengths to RGB colors
function cmap = wavelengthToRGB(lambda)
    % Mapping of wavelengths to RGB colors
    R = zeros(size(lambda));
    G = zeros(size(lambda));
    B = zeros(size(lambda));
    
    % Mapping to the colors of the visible spectrum (simple approximation)
    for i = 1:length(lambda)
        if lambda(i) >= 380 && lambda(i) < 440
            R(i) = -(lambda(i) - 440) / (440 - 380);
            G(i) = 0;
            B(i) = 1;
        elseif lambda(i) >= 440 && lambda(i) < 490
            R(i) = 0;
            G(i) = (lambda(i) - 440) / (490 - 440);
            B(i) = 1;
        elseif lambda(i) >= 490 && lambda(i) < 550
            R(i) = 0;
            G(i) = 1;
            B(i) = -(lambda(i) - 550) / (550 - 490);
        elseif lambda(i) >= 550 && lambda(i) < 590
            R(i) = (lambda(i) - 550) / (590 - 550);
            G(i) = 1;
            B(i) = 0;
        elseif lambda(i) >= 590 && lambda(i) < 650
            R(i) = 1;
            G(i) = -(lambda(i) - 650) / (650 - 590);
            B(i) = 0;
        elseif lambda(i) >= 650 && lambda(i) <= 780
            R(i) = 1;
            G(i) = 0;
            B(i) = 0;
        elseif lambda(i) >= 780 && lambda(i) <= 1000
            % Map from red (R=1) to black (R=0) above 780nm with darker gradient
            R(i) = max(0, 1 - (lambda(i) - 780) / (1000 - 780))^2; % Exponential adjustment for sharper darkening
            G(i) = 0;
            B(i) = 0;
        end
    end
    
    % Normalize RGB values to the range [0, 1]
    R = max(0, min(1, R));
    G = max(0, min(1, G));
    B = max(0, min(1, B));
    
    % Combine the RGB components into a color matrix
    cmap = [R', G', B'];
end