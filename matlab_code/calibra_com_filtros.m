clc; clear all; close all;
global WaveLenght matrix_GSCM;

%% Configuração inicial

% Carrega matriz de correção espectral (Golden Spec Correction Matrix)
WaveLenght = 380:1000;
matrix_GSCM = xlsread("General_Spec_Corr_Matrix.xlsx");
matrix_GSCM = matrix_GSCM(2:end, 2:end);

% Define experimento
experiment = "RGBW_100";

% Carrega dados sem difusor (referência)
load("calibration_tests/sensor_measurements/" + experiment + "-100_amostras.mat");
sensor_values_sem_diff = sum(buffer_sensor_values,1)./99;

% Carrega dados com difusor (sensor a ser calibrado)
load("calibration_tests/sensor_measurements/" + experiment + "_difo1-100_amostras.mat");
sensor_values_com_diff = sum(buffer_sensor_values,1)./99;

% Aplica correção de ganho
gain_correction_256 = 0.987308373187068;
gain_correction_512 = 0.959349243600411;
gain_correction = [gain_correction_512(ones(1,6)), gain_correction_256(ones(1,4))];

sensor_values_sem_diff = sensor_values_sem_diff .* gain_correction;
sensor_values_com_diff = sensor_values_com_diff .* gain_correction;

% Gera SPD de referência (sem difusor)
SPD_ref = matrix_GSCM * sensor_values_sem_diff';

%% Otimização com filtros FIR

% Parâmetros
n_channels = 10;
N = 15; % Ordem do filtro
initial_guess = ones(n_channels * N, 1); % Chute inicial

% Define função objetivo para lsqnonlin
options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt', ...
    'MaxIterations', 500, 'FunctionTolerance', 1e-9);

objective = @(h) fir_filter_objective(h, sensor_values_com_diff, SPD_ref, matrix_GSCM, N);

% Roda otimização
best_h = lsqnonlin(objective, initial_guess, [], [], options);

%% Aplica os filtros encontrados para reconstruir o espectro

sensor_filtered = zeros(n_channels,1);
for k = 1:n_channels
    idx = (k-1)*N + (1:N);
    h_k = best_h(idx)';
    
    % Gera sinal constante prolongado com valor do canal
    x_sim = repmat(sensor_values_com_diff(k), 1, 3*N);
    
    % Aplica filtro com filtfilt para evitar deslocamento
    y_filt = filtfilt(h_k, 1, x_sim);
    
    % Pega valor final como saída estável do canal filtrado
    sensor_filtered(k) = y_filt(end);
end

% Reconstrução do espectro usando os canais filtrados
SPD_filtered = matrix_GSCM * sensor_filtered;

%% Normalização e comparação

SPD_ref = SPD_ref / max(SPD_ref);
SPD_filtered = SPD_filtered / max(SPD_filtered);
SPD_original = matrix_GSCM * sensor_values_com_diff';
SPD_original = SPD_original / max(SPD_original);

% Plot para comparação
figure;
hold on;
plot(WaveLenght, SPD_ref, '--k', 'LineWidth', 1.5);
plot(WaveLenght, SPD_original, 'r', 'LineWidth', 1.2);
plot(WaveLenght, SPD_filtered, 'b', 'LineWidth', 1.2);
legend('Sem difusor (referência)', 'Com difusor (original)', 'Com difusor + FIR calibrado');
xlabel("Wavelength [nm]");
ylabel("Relative Magnitude");
title("Calibração com Filtros FIR (N = " + N + ")");
grid on;

%% Função objetivo
function err = fir_filter_objective(h, sensor_raw, SPD_target, matrix_GSCM, N)
    n_channels = length(sensor_raw);
    sensor_filtered = zeros(n_channels,1);
    
    for k = 1:n_channels
        idx = (k-1)*N + (1:N);
        h_k = h(idx)';

        % Sinal constante artificial para simular o canal
        x_sim = repmat(sensor_raw(k), 1, 3*N);
        y_filt = filtfilt(h_k, 1, x_sim);

        sensor_filtered(k) = y_filt(end); % Valor estabilizado
    end

    SPD_est = matrix_GSCM * sensor_filtered;
    SPD_est = SPD_est / max(SPD_est);
    SPD_target = SPD_target / max(SPD_target);

    err = SPD_est - SPD_target;  % Erro espectral
end
