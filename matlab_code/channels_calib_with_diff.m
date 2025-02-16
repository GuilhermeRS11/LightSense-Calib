
%% Inicialização e Configuração de Ambiente
clc;        % Limpa a janela de comando
clear all;  % Limpa todas as variáveis
close all;  % Fecha todas as figuras

% Inicialização dos valores de calibração baseados em dados padrão da OSRAM
F1 = 0.011057; F2 = 0.019664; F3 = 0.028302; F4 = 0.032492;
F5 = 0.034778; F6 = 0.034016; F7 = 0.039658; F8 = 0.043168;
Clear = 0.127734; NIR = 0.051806;

% Armazenamento global de variáveis para uso em funções
global sensor_values matrix_GSCM;
sensor_values = [F1, F2, F3, F4, F5, F6, F7, F8, NIR];

% Configuração dos Comprimentos de Onda e Matriz de Correção

WaveLenght = 380:1:1000; % Comprimentos de onda total disponíveis

% Carregamento e configuração da matriz de correção espectral
matrix_GSCM = xlsread("Tabelas/General_Spec_Corr_Matrix.xlsx");
matrix_GSCM = matrix_GSCM(2:end, 2:end);

% Aplicação de correções de ganho para ajuste dos valores dos sensores
gain_correction_256 = 0.987308373187068; % Ganho para 256
gain_correction_512 = 0.959349243600411; % Ganho para 512
gain_correction = [gain_correction_512(ones(1,6)), gain_correction_256(ones(1,4))];

% Carrega os coeficientes medios definidos na calibração dos canais
load('calibr_coeficientes.mat');

%% Ajuste dos Contadores Básicos Utilizando o Espectro Reconstruído

% Define os ensaios que serão analisados
ensaios_names = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", ...
                 "GW_100", "RGB_100", "RGBW_100"];

% Configuração para o número de objetivos (tipos de erro) e difusores
n_func_objetivo = 4; % EQM, EQMP, EAM, REQM
n_difo = 2;          % Difusor 1 e Difusor 2
n_canais_sensor = 10; % Número de canais do sensor

% Preparação para armazenar os coeficientes médios obtidos pela calibração
mean_coeficientes = zeros(n_func_objetivo, n_difo, n_canais_sensor);

% Configurações das opções de otimização para o algoritmo lsqnonlin
options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt',... 
    'OptimalityTolerance', 1e-15, 'FunctionTolerance', 1e-50, 'MaxIterations', 1000, ...
    'StepTolerance', 1e-15, 'MaxFunctionEvaluations', 2.000000e+04);

% Iteração sobre todos os ensaios
for i = 1:length(ensaios_names)
    for j = 1:n_difo
        % Carrega os dados do sensor com difusor
        load(append("Ensaios_completos/Sensor_Matlab/", ensaios_names(i), sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_Difo = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
        
        % Carrega os dados sem difusor para comparação
        load(append("Ensaios_completos/Sensor_Matlab/", ensaios_names(i), "-100_amostras.mat"));
        sensor_values_semDifo = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
        
        % Dados originais obtidos da combinação de matriz de correção e valores do sensor sem difusor
        dados_originais = matrix_GSCM * sensor_values_semDifo';
        
        % Função que reconstrói o espectro a partir dos parâmetros ajustados
        funcaoSPD = @(params) spectrum_restore(params, sensor_values_Difo);
        
        % Iteração sobre cada função objetivo
        for k = 1:n_func_objetivo
            % Seleção da função objetivo baseada no tipo de erro
            switch k
                case 1 % EQM
                    funcao_objetivo = @(params) mean((funcaoSPD(params) - dados_originais).^2);
                case 2 % EQMP
                    funcao_objetivo = @(params) mean((funcaoSPD(params) - dados_originais).^2 .* dados_originais');
                case 3 % EAM
                    funcao_objetivo = @(params) mean(abs(funcaoSPD(params) - dados_originais));
                case 4 % REQM
                    funcao_objetivo = @(params) sqrt(mean((funcaoSPD(params) - dados_originais).^2));
            end
            
            % Chute inicial e realização do ajuste dos parâmetros
            parametros_iniciais = ones(1, n_canais_sensor);
            parametros_ajustados = lsqnonlin(funcao_objetivo, parametros_iniciais, [], [], options);
            
            % Armazenamento dos coeficientes ajustados
            mean_coeficientes(k, j, :) = parametros_ajustados;
        end
    end
end

%% Cálculo do R² para cada ensaio usando cada um dos métodos

% Definição dos ensaios e configurações dos testes
ensaios_names = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", "GW_100", "RGB_100", "RGBW_100"];
ensaios_labels = ["R", "G", "B", "W", "RB", "RG", "GW", "RGB", "RGBW"];
numEnsaios = numel(ensaios_names);
n_func_objetivo = 5;  % EQM, EQMP, EAM, REQM, sem calibração
n_difusores = 2;      % Difusor 1 e Difusor 2
R2_results = zeros(numEnsaios, n_func_objetivo, n_difusores); % Matriz para armazenar os resultados de R²

% Loop para processar cada ensaio
for i = 1:numEnsaios
    ensaio = ensaios_names(i);
    
    % Carregar dados do sensor sem difusor para comparação
    load(append("Ensaios_completos/Sensor_Matlab/", ensaio, "-100_amostras.mat"));
    sensor_values_semDifo = sum(buffer_sensor_values, 1) ./ 99 .* gain_correction .* calibr_coeficientes;
    
    % Obtenção do espectro reconstruído sem uso de difusor para referência
    reconstructed_semDifo = spectrum_restore(ones(1,10), sensor_values_semDifo);
    y_mean = mean(reconstructed_semDifo);
    SS_tot = sum((reconstructed_semDifo - y_mean).^2);

    % Processamento para cada configuração de difusor
    for j = 1:n_difusores
        load(append("Ensaios_completos/Sensor_Matlab/", ensaio, sprintf("_difo%d-100_amostras.mat", j)));
        sensor_values_Difo = sum(buffer_sensor_values,1)./99 .* gain_correction .* calibr_coeficientes;

        % Cálculo de R² para o caso sem ajuste de calibração
        reconstructed = spectrum_restore(ones(1,10), sensor_values_Difo);
        SS_res = sum((reconstructed_semDifo - reconstructed).^2);
        R2_results(i, 1, j) = 1 - SS_res / SS_tot;

        % Cálculo de R² para cada método de calibração
        for k = 2:n_func_objetivo
            mean_coef = squeeze(mean_coeficientes(k-1, j, :))'; % Ajuste para alinhar com as matrizes
            reconstructed = spectrum_restore(mean_coef, sensor_values_Difo);
            SS_res = sum((reconstructed_semDifo - reconstructed).^2);
            R2_results(i, k, j) = 1 - SS_res / SS_tot;
        end
    end
end

%% Função que faz a recontrução do espectro através da matriz do Golden Device
function SPD_reconstruida = spectrum_restore(coeficientes_calibracao, sensor_values)
    % Importa variáveis globais necessárias para a função
    global matrix_GSCM;
    
    % Definição dos limites para os coeficientes de calibração
    limites_inferiores = ones(1, 10) * -1;
    limites_superiores = ones(1, 10) * 10;
    
    % Verificação se os coeficientes de calibração estão dentro dos limites estabelecidos
    if any(coeficientes_calibracao < limites_inferiores) || any(coeficientes_calibracao > limites_superiores)
        disp("Coeficiente informado fora dos limites permitidos");
        % Retorna um vetor de erro se os coeficientes estiverem fora dos limites
        SPD_reconstruida_norm = ones(size(WaveLenght)) * 1e10;
        return;
    end

    % Se os coeficientes estão dentro dos limites, realiza o cálculo
    % Multiplica os valores do sensor pelos coeficientes de calibração
    sensor_values_calibr = sensor_values .* coeficientes_calibracao;
    % Multiplica a matriz de correção espectral pelos valores do sensor calibrados
    SPD_reconstruida = matrix_GSCM * sensor_values_calibr';
    % Normaliza o espectro reconstruído pelo seu valor máximo
    SPD_reconstruida_norm = SPD_reconstruida / max(SPD_reconstruida);
end
