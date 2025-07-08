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

%% Frequency response analysis of the GSCM matrix rows
% Each row of the GSCM matrix is interpreted as a FIR filter with 10 coefficients (1 per channel)

[num_lines, num_channels] = size(matrix_GSCM);

f = figure;
hold on;

% Plot frequency response of all rows
for i = 1:num_lines
    h = matrix_GSCM(i, :);  % Coefficients of row i
    [H, w] = freqz(h, 1, 512);  % Frequency response with 512 points
    plot(w/pi, abs(H));  % Normalized frequency (0 to 1)
end

xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');
title('Frequency responses of the GSCM matrix rows');
grid on;
exportgraphics(f, 'saved_images/GSCM_line_freq_response.png', 'ContentType', 'vector');

%% Frequency analysis of the spectral error
% Evaluates the frequency composition of the error (difference between reference and reconstructed SPD)

% Compute pointwise error
erro_spectrum = sphere_SPD' - SPD_reconstructed_adjusted;

% Apply FFT to the error
erro_fft = fft(erro_spectrum);
erro_fft_mag = abs(erro_fft);
erro_fft_mag = erro_fft_mag(1:floor(end/2)); % keep only the positive half (real spectrum)
freqs = linspace(0, 0.5, length(erro_fft_mag)); % Normalized frequency (0 to 0.5)

% Plot frequency spectrum of the error
f = figure;
plot(freqs, erro_fft_mag, 'k', 'LineWidth', 1.5);
xlabel('Normalized Frequency');
ylabel('|FFT(error)|');
title('Frequency analysis of the spectral error');
exportgraphics(f, 'saved_images/Spectral_frequency_error_analisys.png', 'ContentType', 'vector');
grid on;

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

% Gráfico de barras com lsqnonlin e EQM (MSE) apenas
barWidth = 0.25;

barColors = [
    [0.3 0.3 0.3];  % Sem calibração
    [1 0.3 0.3];    % EQM (MSE)
];

f = figure;
hold on;

positions = 1:num_experiments;

bar(positions - barWidth/2, values_no_calib, barWidth, 'FaceColor',  barColors(1, :));  % Sem calibração
bar(positions + barWidth/2, R2_MSE_values(1,:), barWidth, 'FaceColor', barColors(2, :)); % EQM (lsqnonlin)

set(gca, 'XTick', positions, 'XTickLabel', experiment_labels);
xtickangle(45);
legend({'Sem calibração', 'EQM - lsqnonlin'}, 'Location', 'bestoutside');
title("EQM - lsqnonlin");
xlabel('Experimento');
ylabel('R²');
grid on;
axis tight;
min_value = min([min(values_no_calib), min(R2_MSE_values(1,:))]);
ylim([min_value - abs(min_value * 0.2), 1.1]);
hold off;



%set(gcf, 'Position', get(0, 'Screensize'));  % Maximizes the figure window
exportgraphics(f,'saved_images/All_R2_EQM_lsqninlin.png','ContentType','vector');

%% Reconstructed spectrum using the average coefficient of each method

% Define the experiment to be analyzed. Within the standard of "Complete Experiments"
experiment = "W_100";
Nome_grafico = "White";
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
plot(WaveLenght, sphere_SPD,'--k','linew',2); % Referência da esfera
title(Nome_grafico);
hold on;

plot(WaveLenght, reconstruction_no_diff,'Color', barColors(1,:),'lineWidth', 1); % Sem calibração
plot(WaveLenght, SPD_MSE,'Color', barColors(2,:),'lineWidth', 1); % EQM - lsqnonlin

minValue = min([min(reconstruction_no_diff), min(SPD_MSE), min(sphere_SPD)]);
maxValue = max([max(reconstruction_no_diff), max(SPD_MSE), max(sphere_SPD)]);

for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*(0.15), minValue - maxValue*(0.15), minValue - maxValue*(0.05), minValue - maxValue*(0.05)];
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

xlabel("Comprimento de onda [nm]");
ylabel("Magnitude relativa");
grid on;
legend(["Referência (esfera)", "Sem calibração", "EQM - lsqnonlin"], 'location', 'northeast');
ylim([minValue - maxValue*(0.15), maxValue*(1.1)]);
xlim([380, 1000]);


exportgraphics(f,append('saved_images/SPD_reconstructed_methods_',experiment,'.png'),'ContentType','vector');


%% Post-processing of the reconstructed SPD (spectral smoothing)
% Applies smoothing filters to reduce undesired fluctuations

% Apply Savitzky-Golay filter (order 2, window size 11)
SPD_smooth_sgolay = sgolayfilt(SPD_MSE, 2, 11);

% Apply moving average filter as an alternative
SPD_smooth_movmean = movmean(SPD_MSE, 11);

% Compute RMSE before and after smoothing
rmse_original = sqrt(mean((SPD_MSE - sphere_SPD').^2));
rmse_sgolay = sqrt(mean((SPD_smooth_sgolay - sphere_SPD').^2));
rmse_movmean = sqrt(mean((SPD_smooth_movmean - sphere_SPD').^2));

disp("RMSE - Adjusted (no filter): " + rmse_original);
disp("RMSE - Smoothed (Savitzky-Golay): " + rmse_sgolay);
disp("RMSE - Smoothed (Moving Average): " + rmse_movmean);

% Plotting to compare results
f = figure;
plot(WaveLenght, sphere_SPD, '--k', 'LineWidth', 2);
hold on;
plot(WaveLenght, SPD_MSE, 'b', 'LineWidth', 1);
plot(WaveLenght, SPD_smooth_sgolay, 'r', 'LineWidth', 1);
plot(WaveLenght, SPD_smooth_movmean, 'g', 'LineWidth', 1);
legend('Reference (integrating sphere)', 'Adjusted reconstruction', ...
       'Savitzky-Golay', 'Moving average', 'Location', 'northeast');
xlabel("Wavelength [nm]");
ylabel("Relative Magnitude");
grid on;
title("Spectral smoothing comparison of the reconstructed SPD");
exportgraphics(f, 'saved_images/SPD_smoothing_comparison_RGB_100.png', 'ContentType', 'vector');

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