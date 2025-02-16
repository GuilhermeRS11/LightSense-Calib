

%% Inicialização e Configuração de Ambiente
clc;        % Limpa a janela de comando
clear all;  % Limpa todas as variáveis
close all;  % Fecha todas as figuras

% Inicialização dos valores de calibração baseados em dados padrão da OSRAM
F1 = 0.011057; F2 = 0.019664; F3 = 0.028302; F4 = 0.032492;
F5 = 0.034778; F6 = 0.034016; F7 = 0.039658; F8 = 0.043168;
Clear = 0.127734; NIR = 0.051806;

% Armazenamento global de variáveis para uso em funções
global sensor_values WaveLenght matrix_GSCM;
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

% Definindo os ensaios e as configurações iniciais
ensaios_names = ["R_100", "G_100", "B_100", "W_100", "RB_100", "RG_100", "GW_100", "RGB_100", "RGBW_100"];
ensaios_labels = ["R", "G", "B", "W", "RB", "RG", "GW", "RGB", "RGBW"];
numEnsaios = numel(ensaios_names);
n_func_otimizadoras = 2;  % lsqnonlin (1), particleswarm (2)
n_func_objetivo = 4;       % EQM, EQMP, EAM, REQM
n_canais_sensor = 10;

%% Método 1 (Apenas 1 ensaio) - Ajustar basic_counts através do espectro recontruido
% Os fragmentos de código comentados podem ser utilizados para analisar
% outras funções objetivo e outras funções de otimização

ensaio = 'RB_100';
global SPD_esfera;

% Carrega os dados dos ensaios do sensor e da esfera integradora
load(append("Ensaios_completos/Sensor_Matlab/", ensaio,"-100_amostras.mat"));
sensor_values_semDifo = sum(buffer_sensor_values,1)./99;

% Aplica a calibração de ganho
sensor_values_semDifo = sensor_values_semDifo .* gain_correction;

% Adapta os valores da esfera de 800 a 1100 para serem zero, já que não temos informção deles
SPD_esfera = zeros(1, length(WaveLenght)); 
load(append("Ensaios_completos/Esfera/Espectro/", ensaio,"_S.mat"));
SPD_esfera(1:421) = data.Wrad_relative;
coef_correcao_semDifo = zeros(length(sensor_values_semDifo));

% Função que gera a SPD reconstruída a partir dos parâmetros
funcaoSPD = @(params) spectrum_restore(params, sensor_values_semDifo);

% Funções objetivo para minimizar a diferença entre a SPD reconstruída e a SPD original
%EQM (Erro Quadrático Médio) não ponderado
funcao_objetivo = @(params) mean((funcaoSPD(params) - SPD_esfera').^2); 

% EQM (Erro Quadrático Médio) ponderado
%funcao_objetivo = @(params) mean((funcaoSPD(params) - SPD_esfera').^2 .* SPD_esfera');

% EAM (Erro Absoluto Médio) 
%funcao_objetivo = @(params) mean(abs(funcaoSPD(params) - SPD_esfera'));

% REQM (Raiz do Erro Quadrático Médio)
% funcao_objetivo = @(params) sqrt(mean((funcaoSPD(params) - SPD_esfera').^2));

% Chute inicial para os 10 parâmetros ajustáveis
parametros_iniciais = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

%Realize o ajuste dos parâmetros usando lsqnonlin
options = optimoptions('lsqnonlin','Display','iter'); % Para exibir informações sobre o processo de otimização
options.Algorithm = 'levenberg-marquardt';
options.OptimalityTolerance = 1e-15;
options.FunctionTolerance = 1e-50;
options.MaxIterations = 1000; % Defina o número máximo de iterações
options.StepTolerance = 1e-15; % Defina a tolerância do passo
options.MaxFunctionEvaluations = 1e4;
parametros_ajustados = lsqnonlin(funcao_objetivo, parametros_iniciais, [], [], options);

% lb = 0.01 * ones(1, 10);  % Limite inferior para todos os coeficientes
% ub = 10 * ones(1, 10);   % Limite superior para todos os coeficientes
% options = optimoptions('particleswarm', 'SwarmSize', 30, 'HybridFcn', @fmincon, 'Display', 'iter');
% [parametros_ajustados, fval] = particleswarm(funcao_objetivo, 10, lb, ub, options);

SPD_recontruida_original = funcaoSPD(parametros_iniciais);
% Utiliza os parâmetros ajustados para obter a SPD reconstruída otimizada
SPD_reconstruida_ajustada = funcaoSPD(parametros_ajustados);
coef_correcao_semDifo = parametros_ajustados;      

EQM = mean((SPD_recontruida_original - SPD_esfera').^2); % Erro quadrado médio
EQMP = mean((SPD_recontruida_original - SPD_esfera').^2.*SPD_esfera'); % Erro quadrado médio ponderado
disp(append("EQM no inicio = ",num2str(EQM)));
disp(append("EQMP no inicio = ",num2str(EQMP)));

EQM = mean((SPD_reconstruida_ajustada - SPD_esfera').^2); % Erro quadrado médio
EQMP = mean((SPD_reconstruida_ajustada - SPD_esfera').^2.*SPD_esfera'); % Erro quadrado médio ponderado
disp(append("EQM no fim = ",num2str(EQM)));
disp(append("EQMP no fim = ",num2str(EQMP)));

% Plote dos resultados para avaliação
cmap = wavelengthToRGB(WaveLenght);

f = figure;
h1 = plot(WaveLenght, SPD_esfera,'--k','linew',2);
%title("Calibração dos canais");
hold on;
h2 = plot(WaveLenght, SPD_recontruida_original, 'r', 'LineWidth', 1);
h3 = plot(WaveLenght, SPD_reconstruida_ajustada, 'b', 'LineWidth', 1);

minValue = min(min(SPD_recontruida_original), min(SPD_reconstruida_ajustada));
maxValue = 1;
% Adicionar retângulo colorido para representar as cores do espectro
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*(0.15), minValue - maxValue*(0.15), minValue - maxValue*(0.05), minValue - maxValue*(0.05)]; % Altura do retângulo (ajuste conforme necessário)
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

xlabel("Comprimento de onda [nm]");
ylabel("Magnitude relativa");
grid on; 

leg = legend([h1, h2, h3], ["SPD da esfera", "Reconstrução", "Reconstrução ajustada"], 'Location', 'northeast');
ylim([minValue - maxValue*(0.15), maxValue*(1.1)]);
xlim([WaveLenght(1) WaveLenght(end)]);

exportgraphics(f,append('Saved_images/CalibraCoeficientesIdeais_',ensaio,'.pdf'),'ContentType','vector');


%% Método 1 (Todos os ensaios) - Ajustar basic_counts através do espectro recontruido

mean_coeficientes = zeros(n_func_otimizadoras, n_func_objetivo, n_canais_sensor);
std_coeficientes = zeros(n_func_otimizadoras, n_func_objetivo, n_canais_sensor);
coeficientes_correcao = zeros(n_func_otimizadoras, n_func_objetivo, numEnsaios, n_canais_sensor);
tic;

% Processa todos os dados para o numero de funções otimizadoras
for k = 1:n_func_otimizadoras
    % Processa todos os dados para o numero de funções objetivo
    for j = 1:n_func_objetivo
        % Faz a calibração de todos os ensaios solicitados
        for i = 1:numEnsaios  

            % Carrega os dados dos ensaios do sensor e da esfera integradora
            load(append("Ensaios_completos/Sensor_Matlab/", ensaios_names(i),"-100_amostras.mat"));
            sensor_values_semDifo = sum(buffer_sensor_values,1)./99;

            % Aplica a calibração de ganho
            sensor_values_semDifo = sensor_values_semDifo .* gain_correction;

            % Adapta os valores da esfera de 800 a 1100 para serem zero, já que não temos informção deles
            SPD_esfera = zeros(1, length(WaveLenght)); 
            load(append("Ensaios_completos/Esfera/Espectro/", ensaios_names(i),"_S.mat"));
            SPD_esfera(1:421) = data.Wrad_relative;

            % Função que gera a SPD reconstruída a partir dos parâmetros
            funcaoSPD = @(params) spectrum_restore(params, sensor_values_semDifo);
                     
            % Função objetivo para minimizar a diferença entre a SPD reconstruída e a SPD original
            if j == 1
                % EQM (Erro Quadrático Médio) não ponderado
                funcao_objetivo = @(params) mean((funcaoSPD(params) - SPD_esfera').^2);
            
            elseif j == 2
                % EQMP (Erro Quadrático Médio Ponderado) 
                funcao_objetivo = @(params) mean(((funcaoSPD(params) - SPD_esfera').^2 ).* SPD_esfera');

            elseif j == 3
                % EAM (Erro Absoluto Médio) 
                funcao_objetivo = @(params) mean(abs(funcaoSPD(params) - SPD_esfera'));
                
            else    
                % REQM (Raiz do Erro Quadrático Médio)
                funcao_objetivo = @(params) sqrt(mean((funcaoSPD(params) - SPD_esfera').^2));
            end
            
            % Chute inicial para os 10 parâmetros ajustáveis
            parametros_iniciais = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

            if k == 1
                % Realiza a otimização dos parametros com o lsqnonlin
                options = optimoptions('lsqnonlin','Display','iter'); % Para exibir informações sobre o processo de otimização
                options.Algorithm = 'levenberg-marquardt';
                options.OptimalityTolerance = 1e-15;
                options.FunctionTolerance = 1e-50;
                options.MaxIterations = 1000; % Defina o número máximo de iterações
                options.StepTolerance = 1e-15; % Defina a tolerância do passo
                options.MaxFunctionEvaluations = 1e4;
                parametros_ajustados = lsqnonlin(funcao_objetivo, parametros_iniciais, [], [], options);

            else
                % Realiza a otimização dos parametros com o particleswarm
                lb = 0.01 * ones(1, 10);  % Limite inferior para todos os coeficientes
                ub = 10 * ones(1, 10);   % Limite superior para todos os coeficientes
                options = optimoptions('particleswarm', 'SwarmSize', 30, 'HybridFcn', @fmincon, 'Display', 'iter');
                [parametros_ajustados, fval] = particleswarm(funcao_objetivo, 10, lb, ub, options);
            end
            
            coeficientes_correcao(k, j, i, :) = parametros_ajustados;  
        end
        % Calcula a média ao longo da dimensão dos ensaios para cada canal
        mean_coeficientes(k, j, :) = mean(squeeze(coeficientes_correcao(k, j, :, :)), 1);
        std_coeficientes(k, j, :) = std(squeeze(coeficientes_correcao(k, j, :, :)), 1);
    end
end
toc;

%% Calculo do R² para cada ensaio usando cada um dos metodos

% Número de ensaios
numEnsaios = numel(ensaios_names);

% Inicialização correta das variáveis
R2_EQM_values = zeros(n_func_otimizadoras, numEnsaios);
R2_EQMP_values = zeros(n_func_otimizadoras, numEnsaios);
R2_EAM_values = zeros(n_func_otimizadoras, numEnsaios);
R2_REQM_values = zeros(n_func_otimizadoras, numEnsaios);
values_no_calib = zeros(1, numEnsaios);
no_Coeficients = ones(1, 10);

% Loop para calcular as métricas
for j = 1:n_func_otimizadoras
    for i = 1:numEnsaios
        ensaio = ensaios_names(i);
        % Carregar dados do sensor
        load(append("Ensaios_completos/Sensor_Matlab/", ensaio, "-100_amostras.mat"));
        sensor_values_semDifo = sum(buffer_sensor_values,1)./99;
        sensor_values_semDifo = sensor_values_semDifo .* gain_correction;
    
        % Carregar e ajustar dados da esfera
        SPD_esfera = zeros(1, length(WaveLenght));
        load(append("Ensaios_completos/Esfera/Espectro/", ensaio, "_S.mat"));
        SPD_esfera(1:421) = data.Wrad_relative;
        
        % Metrica para avalisar o resempenho dos metodos atraves de R²
        y_mean = mean(SPD_esfera);
        SS_tot = sum((SPD_esfera - y_mean).^2);
        
        % Reconstrução sem calibração
        reconstructed_no_calib = spectrum_restore(no_Coeficients, sensor_values_semDifo);
        SS_res_no_calib = sum((SPD_esfera' - reconstructed_no_calib).^2);
        values_no_calib(1,i) = 1 - SS_res_no_calib / SS_tot;
    
        % Reconstruir usando coeficientes e calcular métricas
        reconstructed = spectrum_restore(squeeze(mean_coeficientes(j, 1, :))', sensor_values_semDifo);
        SS_res = sum((SPD_esfera' - reconstructed).^2);
        R2_EQM_values(j,i) = 1 - SS_res / SS_tot;
        
        reconstructed = spectrum_restore(squeeze(mean_coeficientes(j, 2, :))', sensor_values_semDifo);
        SS_res = sum((SPD_esfera' - reconstructed).^2);
        R2_EQMP_values(j,i) = 1 - SS_res / SS_tot;
         
        reconstructed = spectrum_restore(squeeze(mean_coeficientes(j, 3, :))', sensor_values_semDifo);
        SS_res = sum((SPD_esfera' - reconstructed).^2);
        R2_EAM_values(j,i) = 1 - SS_res / SS_tot;
         
        reconstructed = spectrum_restore(squeeze(mean_coeficientes(j, 4, :))', sensor_values_semDifo);
        SS_res = sum((SPD_esfera' - reconstructed).^2);
        R2_REQM_values(j,i) = 1 - SS_res / SS_tot;
    end
end

% Plotando as barras incluindo espectro sem calibração
barWidth = 0.15;  % Ajuste para acomodar uma barra extra
% Cores especificadas para cada barra
barColors = [
    [0.3 0.3 0.3];  % Cinza escuro para 'Sem calibração'
    [1 0.3 0.3];    % Vermelho suave para EQM
    [0.5 1 0.5];    % Verde suave para EQMP
    [1 0.5 1];      % Magenta suave para EAM
    [0.5 0.5 1];    % Azul suave para REQM
];
% Reconstruir usando coeficientes de todos os métodos para este ensaio
positions = (1:numEnsaios) + (i-1)*barWidth*4 - barWidth*(numEnsaios/2-0.5); % Posições ajustadas para as barras
R2_values_lsqnonlin = [R2_EQM_values(1,:) R2_EQMP_values(1,:) R2_EAM_values(1,:) R2_REQM_values(1,:)];
R2_values_particleswarm = [R2_EQM_values(2,:) R2_EQMP_values(2,:) R2_EAM_values(2,:) R2_REQM_values(2,:)];

f = figure;
for j = 1:n_func_otimizadoras
    subplot(n_func_otimizadoras,1,j);
    hold on;
    bar(positions - 2 * barWidth, values_no_calib, barWidth, 'FaceColor',  barColors(1, :));  % Preto para não calibrado
    bar(positions - 1 * barWidth, R2_EQM_values(j,:), barWidth, 'FaceColor', barColors(2, :));
    bar(positions, R2_EQMP_values(j,:), barWidth, 'FaceColor', barColors(3, :));
    bar(positions + 1 * barWidth, R2_EAM_values(j,:), barWidth, 'FaceColor', barColors(4, :));
    bar(positions + 2 * barWidth, R2_REQM_values(j,:), barWidth, 'FaceColor', barColors(5, :));
    
    % Configurações do gráfico
    set(gca, 'XTick', positions, 'XTickLabel', ensaios_labels);
    xtickangle(45);
    legend({'Sem calibração', 'EQM', 'EQMP', 'EAM', 'REQM'}, 'Location', 'bestoutside');
    title([titles_otimizadora(j)]);
    xlabel('Ensaio');
    ylabel('R²');
    grid on;
    axis tight;
    min_value = min(min(R2_values_lsqnonlin), min(R2_values_particleswarm));
    ylim([min_value - abs(min_value * 0.5), 1.1]);
    hold off;
end

% Exportar os coeficientes medios
calibr_coeficientes = squeeze(mean_coeficientes(1, 1, :))'; % Exporta os EQM que obitiveram melhor resultado
save('calibr_coeficientes.mat', 'calibr_coeficientes');

set(gcf, 'Position', get(0, 'Screensize'));  % Maximiza a janela da figura
exportgraphics(f,'Saved_images/R2_todas.pdf','ContentType','vector');

%% Espectro reconstruido utilizando o coeficiente médio de cada método

% Define o ensaio que será analisado. Dentro do padrão dos "Ensaios Completos"
ensaio = "RGBW_100";
Nome_grafico = "RGBW";
load(append("Ensaios_completos/Sensor_Matlab/",ensaio,"-100_amostras.mat"));
file = append("Ensaios_completos/Esfera/Espectro/",ensaio,"_S.txt")';

load(append("Ensaios_completos/Sensor_Matlab/",ensaio,"-100_amostras.mat"));
sensor_values_semDifo = sum(buffer_sensor_values,1)./99;
sensor_values_semDifo = sensor_values_semDifo .* gain_correction;

% Adapta os valores da esfera de 800 a 1100 para serem zero, já que não temos informção deles
SPD_esfera = zeros(1, length(WaveLenght)); 
load(append("Ensaios_completos/Esfera/Espectro/", ensaio,"_S.mat"));
SPD_esfera(1:421) = data.Wrad_relative;


reconstrucao_semDifo = matrix_GSCM * sensor_values_semDifo';
reconstrucao_semDifo = reconstrucao_semDifo ./ max(reconstrucao_semDifo);

SPD_EQM = spectrum_restore(squeeze(mean_coeficientes(1, 1, :))', sensor_values_semDifo);
SPD_EQMP = spectrum_restore(squeeze(mean_coeficientes(1, 2, :))', sensor_values_semDifo);
SPD_EAM = spectrum_restore(squeeze(mean_coeficientes(1, 3, :))', sensor_values_semDifo);
SPD_REQM = spectrum_restore(squeeze(mean_coeficientes(1, 4, :))', sensor_values_semDifo);

% Para mapear as cores do retangulo que mostra a cor do espectro
cmap = wavelengthToRGB(WaveLenght);

f = figure;
plot(WaveLenght, SPD_esfera,'--k','linew',2);
title(Nome_grafico);
hold on;

plot(WaveLenght, reconstrucao_semDifo,'Color', barColors(1,:),'lineWidth', 1);
plot(WaveLenght, SPD_EQM,'Color', barColors(2,:),'lineWidth', 1);
plot(WaveLenght, SPD_EQMP,'Color', barColors(3,:),'lineWidth', 1);
plot(WaveLenght, SPD_EAM,'Color', barColors(4,:),'lineWidth', 1);
plot(WaveLenght, SPD_REQM,'Color', barColors(5,:),'lineWidth', 1);


minValue = min([min(reconstrucao_semDifo), min(SPD_EQM), min(SPD_EQMP), min(SPD_EAM), min(SPD_REQM), min(SPD_esfera)]);
maxValue = max([max(reconstrucao_semDifo), max(SPD_EQM), max(SPD_EQMP), max(SPD_EAM), max(SPD_REQM), max(SPD_esfera)]);
% Adicionar retângulo colorido para representar as cores do espectro
for i = 1:length(WaveLenght)-1
    x = [WaveLenght(i), WaveLenght(i+1), WaveLenght(i+1), WaveLenght(i)];
    y = [minValue - maxValue*(0.15), minValue - maxValue*(0.15), minValue - maxValue*(0.05), minValue - maxValue*(0.05)]; % Altura do retângulo (ajuste conforme necessário)
    fill(x, y, cmap(i,:), 'EdgeColor', 'none');
end

xlabel("Comprimento de onda [nm]");
ylabel("Magnitude relativa");
grid on; 
legend(["SPD da esfera" "Sem calibrar" "EQM" "EQMP" "EAM" "REQM"],'location','northeast');
ylim([minValue - maxValue*(0.15), maxValue*(1.1)]);
xlim([380, 1000]);

exportgraphics(f,append('Saved_images/SPD_recontruido_metodos',ensaio,'.pdf'),'ContentType','vector');

%% Função que faz a recontrução do espectro através da matriz do Golden Device
function SPD_reconstruida_norm = spectrum_restore(coeficientes_calibracao, sensor_values)
    global WaveLenght matrix_GSCM;
    
    % Limites para os coeficientes de calibração
    limites_inferiores = ones(1, 10) * 0.01;
    limites_superiores = ones(1, 10) * 10;
    
    % Checa se os coeficientes estão dentro dos limites
    if any(coeficientes_calibracao < limites_inferiores) || any(coeficientes_calibracao > limites_superiores)
        disp("Coeficiente informado fora dos limites permitidos");
        SPD_reconstruida_norm = ones(length(WaveLenght),1) * 1e10;  % Saída de erro
        return;
    end

    % Se os coeficientes estão dentro dos limites, procede com o cálculo
    sensor_values_calibr = sensor_values .* coeficientes_calibracao;
    SPD_reconstruida = matrix_GSCM * sensor_values_calibr';
    SPD_reconstruida_norm = SPD_reconstruida / max(SPD_reconstruida);  % Normaliza o resultado
      
end

% Função para mapear comprimentos de onda para cores RGB
function cmap = wavelengthToRGB(lambda)
    % Mapeamento de comprimentos de onda para cores RGB
    R = zeros(size(lambda));
    G = zeros(size(lambda));
    B = zeros(size(lambda));
    
    % Mapeamento para as cores do espectro visível (aproximação simples)
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
            % Mapear de vermelho (R=1) para preto (R=0) acima de 780nm com gradiente mais escuro
            R(i) = max(0, 1 - (lambda(i) - 780) / (1000 - 780))^2; % Ajuste exponencial para escurecimento mais acentuado
            G(i) = 0;
            B(i) = 0;
        end
    end
    
    % Normalizar os valores RGB para o intervalo [0, 1]
    R = max(0, min(1, R));
    G = max(0, min(1, G));
    B = max(0, min(1, B));
    
    % Combinar os componentes RGB em uma matriz de cores
    cmap = [R', G', B'];
end