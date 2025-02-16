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
        load(append("calibration_tests/sensor_measurements/", experiment_name(i), sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_with_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;
        
        % Load data without diffuser for comparison
        load(append("calibration_tests/sensor_measurements/", experiment_name(i), "-100_amostras.mat"));
        sensor_values_without_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;
        
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
experiment_names = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", "GW_100", "RGB_100", "RGBW_100"];
experiment_labels = ["R", "G", "B", "W", "RB", "RG", "GW", "RGB", "RGBW"];
num_experiments = numel(experiment_names);
n_opt_func = 5;       % MSE, WMSE, MAE, RMSE, without calibration
n_diffs = 2;      % Diffuser 1 and Diffuser 2
R2_results = zeros(num_experiments, n_opt_func, n_diffs); % Matrix to store R² results

% Loop to process each trial
for i = 1:num_experiments
    experiment = experiment_names(i);
    
    % Load sensor data without diffuser for comparison
    load(append("calibration_tests/sensor_measurements/", experiment, "-100_amostras.mat"));
    sensor_values_without_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;
    
    % Obtain the reconstructed spectrum without using a diffuser for reference
    reconstructed_semDifo = spectrum_restore(ones(1,10), sensor_values_without_diffuser);
    y_mean = mean(reconstructed_semDifo);
    SS_tot = sum((reconstructed_semDifo - y_mean).^2);

    % Processing for each diffuser configuration
    for j = 1:n_diffs
        load(append("calibration_tests/sensor_measurements/", experiment, sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_with_diffuser = sum(buffer_sensor_values,1)./99 .* gain_correction .* calibr_coeficients;

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

% Colors for the bars
barColors = [
    [0.3 0.3 0.3];  % Dark gray for 'No calibration'
    [1 0.3 0.3];    % Soft red for MSE
    [0.5 1 0.5];    % Soft green for WMSE
    [1 0.5 1];      % Soft magenta for MAE
    [0.5 0.5 1];    % Soft blue for RMSE
];

% Adjusting the bar width and group width
barWidth = 0.15;
groupWidth = n_opt_func * barWidth;

% Define font sizes
fontSizeTitle = 18-6;
fontSizeLabels = 18-6;
fontSizeAxis = 18-6;
fontSizeLegend = 18-6;

% Plot for each diffuser
f = figure;

for j = 1:n_diffs
    subplot(n_diffs, 1, j);
    hold on;
    title(sprintf('Difusor %d', j), 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    
    for i = 1:num_experiments
        basePosition = (i - 1) * (groupWidth + barWidth) + barWidth;

        for k = 1:n_opt_func
            barPosition = basePosition + (k - 1) * barWidth;
            bar(barPosition, R2_results(i, k, j), barWidth, 'FaceColor', barColors(k, :));
        end
    end

    set(gca, 'XTick', barWidth*3 + (0:num_experiments-1) * (groupWidth + barWidth), 'XTickLabels', experiment_labels, 'FontSize', fontSizeAxis);
    xtickangle(45);
    ylabel('R²', 'FontSize', fontSizeLabels);
    xlabel('Experiment', 'FontSize', fontSizeLabels);
    legend({'No calibration', 'MSE', 'WMSE', 'MAE', 'RMSE'}, 'Location', 'bestoutside', 'FontSize', fontSizeLegend);
    grid on;
    axis tight;
    min_value = min(min(min(R2_results)));
    ylim([min_value - abs(min_value * 0.5), 1.1]);
    hold off;
end

screenSize = get(0, 'Screensize');
set(gcf, 'Position', [screenSize(1), screenSize(2), screenSize(3)/2, screenSize(4)]);  % Half the width and full height of the screen
exportgraphics(f,'saved_images/Rsquared_eachMethod_eachTrial_diffusers.pdf','ContentType','vector');

%% Overall R² of each method, aiming to simplify the analysis (average)

% Calibration methods and diffusers
methods = {"No calibration", 'MSE', 'WMSE', 'MAE', 'RMSE'};
diffusers = {'Diffuser 1', 'Diffuser 2'};
n_methods = numel(methods);
n_diffs = numel(diffusers);

% Calculating the average R^2 for each method and diffuser
averages = zeros(n_methods, n_diffs);
for j = 1:n_diffs
    for k = 1:n_methods
        averages(k, j) = mean(R2_results(:, k, j));
    end
end

% Colors for each diffuser
barColors = [
    [0 1 0];  % Green for Diffuser 1
    [0 0 1];  % Blue for Diffuser 2
];

barWidth = 0.4;  % Largura das barras

% Creating the chart
f = figure;
hold on;
% title('Average R² by Calibration Method and Diffuser');

% Iterating over each method to create a group of bars
for k = 1:n_methods
    for j = 1:n_diffs
        % Position of the bar in the group
        barPos = k + (j-1) * barWidth;
        bar(barPos, averages(k, j), barWidth, 'FaceColor', barColors(j, :));
    end
end

% Graph settings
set(gca, 'XTick', 1.2:1:(n_methods + 0.2), 'XTickLabels', methods);
ylabel('R²');
xlabel('Calibration Method');
legend(diffusers, 'Location', 'bestoutside');
grid on;
axis tight;
ylim([0 1.1]);
hold off;

exportgraphics(f,'saved_images/Rsquared_Method_eachMethod.pdf','ContentType','vector');

%% Reconstructed spectrum using the average coefficient of each method

% Initial configuration and data loading
experiment = "RGB_100";
load(append("calibration_tests/sensor_measurements/", experiment, "-100_amostras.mat"));
sensor_values_without_diffuser = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;

load(append("calibration_tests/sensor_measurements/", experiment, "_difo1-100_amostras.mat"));
sensor_values_Diffuser1 = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;

load(append("calibration_tests/sensor_measurements/", experiment, "_difo2-100_amostras.mat"));
sensor_values_Diffuser2 = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficients;

% Spectrum without diffuser as reference
reconstructed_noDiff = matrix_GSCM * sensor_values_without_diffuser';

% Colors for each method as specified
barColors = [
    [1 0.3 0.3];    % Soft red for MSE
    [0.5 1 0.5];    % Soft green for WMSE
    [1 0.5 1];      % Soft magenta for MAE
    [0.5 0.5 1];    % Soft blue for RMSE
    [0.3 0.3 0.3];  % Dark gray for 'No calibration'
];

% Define font sizes
fontSizeTitle = 18-6;
fontSizeLabels = 18-6;
fontSizeAxis = 18-6;
fontSizeLegend = 18-6;

% Generate the color map for the wavelength
cmap = wavelengthToRGB(WaveLenght);

% Names for the legend
LegendNames = {'Without Diffuser', 'No Calibration', 'MSE', 'WMSE', 'MAE', 'RMSE'};

f = figure;
% Adjust the dimensions [left, bottom, width, height] as needed
set(f, 'Position', [100, 100, 1200, 800]); % Larger width and height to avoid flattening

% First subplot
subplot(1,2,1);
hold on;
plotHandle(1) = plot(WaveLenght, reconstructed_noDiff, '--k', 'LineWidth', 2);  % Without Diffuser in dark gray
title("Reconstruction for Diffuser 1", 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
xlabel("Wavelength [nm]", 'FontSize', fontSizeLabels);
ylabel("Irradiance [W/m²]", 'FontSize', fontSizeLabels);

% Plot the reconstruction without calibration
reconstructed_no_calib = spectrum_restore(ones(1,10), sensor_values_Diffuser1);
plotHandle(2) = plot(WaveLenght, reconstructed_no_calib, 'Color', barColors(5,:), 'LineWidth', 1);

% Plot reconstructions for each calibration method
for k = 1:4
    reconstructed = spectrum_restore(squeeze(mean_coeficients(k, 1, :))', sensor_values_Diffuser1);
    plotHandle(k + 2) = plot(WaveLenght, reconstructed, 'Color', barColors(k,:), 'LineWidth', 1);
end

% Add color bar
minValue = min([reconstructed_noDiff; reconstructed_no_calib]);
maxValue = max([reconstructed_noDiff; reconstructed_no_calib]);
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*0.15, minValue - maxValue*0.15, minValue - maxValue*0.05, minValue - maxValue*0.05];
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

colorBarBottom = minValue - maxValue*0.15; % Calculate the bottom position of the colored rectangle
ymin = min(min(reconstructed_noDiff), min(reconstructed_no_calib));
ymax = max(max(reconstructed_noDiff), max(reconstructed_no_calib));

legend(plotHandle, LegendNames, 'FontSize', fontSizeLegend);
grid on;
xlim([380, 1000]);
ylim([colorBarBottom ymax + abs(ymax * 0.2)])
hold off;

% Second subplot
subplot(1,2,2);
hold on;
plotHandle(1) = plot(WaveLenght, reconstructed_noDiff, '--k', 'LineWidth', 2, 'DisplayName', 'Without Diffuser');
title("Reconstruction for Diffuser 2", 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
xlabel("Wavelength [nm]", 'FontSize', fontSizeLabels);
ylabel("Irradiance [W/m²]", 'FontSize', fontSizeLabels);

reconstructed = spectrum_restore(ones(1,10), sensor_values_Diffuser2);
plotHandle(2) = plot(WaveLenght, reconstructed, 'Color', barColors(5,:), 'LineWidth', 1, 'DisplayName', 'No Calibration'); % plot without calibration

for k = 1:4
    reconstructed = spectrum_restore(squeeze(mean_coeficients(k, 2, :))', sensor_values_Diffuser2);
    plotHandle(k + 2) = plot(WaveLenght, reconstructed, 'Color', barColors(k, :), 'LineWidth', 1);
end

% Add color bar
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*0.15, minValue - maxValue*0.15, minValue - maxValue*0.05, minValue - maxValue*0.05];
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

legend(plotHandle, LegendNames, 'FontSize', fontSizeLegend);
grid on;
xlim([380, 1000]);
ylim([colorBarBottom ymax + abs(ymax * 0.2)]);
hold off;

% Adjust the font size of the axis values
ax1 = subplot(1,2,1);
ax1.XAxis.FontSize = fontSizeAxis;
ax1.YAxis.FontSize = fontSizeAxis;

ax2 = subplot(1,2,2);
ax2.XAxis.FontSize = fontSizeAxis;
ax2.YAxis.FontSize = fontSizeAxis;

% Limit the number of legends and adjust the axes to show the graph

screenSize = get(0, 'Screensize');
set(gcf, 'Position', [screenSize(1), screenSize(2), screenSize(3)/2, screenSize(4)]);  % Half the width and full height of the screen
exportgraphics(f,append('saved_images/SPD_reconstructed_methods',experiment,'.pdf'),'ContentType','vector');

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
    
    % Normalize the RGB values to the range [0, 1]
    R = max(0, min(1, R));
    G = max(0, min(1, G));
    B = max(0, min(1, B));
    
    % Combine the RGB components into a color matrix
    cmap = [R', G', B'];
end
