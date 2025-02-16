% This code reads data from the AS7341 sensor sent by the Arduino via serial
% and processes the data, adjusting and displaying the results

%% Process sensor data, display on the graph, and calculate illuminance
clc; 
clear all;
arduino  = serialport("COM5", 9600);

sensor_waveLength = zeros(1,10);
sensor_values = zeros(1,10);

ATIME = 0;
ASTEP = 65534;
GAIN512 = 512; % Para F1 - F6
GAIN256 = 256; % Para F7, F8, Clear

tint_ms = (ATIME + 1) * (ASTEP + 1) * (2.78/1000);

%[0.040000000000000,0.120000000000000,0.100000000000000,0.180000000000000,0.130000000000000,0.100000000000000,0.320000000000000,0.210000000000000,1.190000000000000,0.090000000000000];


% With an offset , the basic counts will be reduced, Offset can be DARK CURRENTS or other 
% constant vlaus those will affects the measurements permantly,
corr_sensor_offset = [0.00196979, 0.00724927, 0.00319381, 0.001314659, 0.001468153, 0.001858105, 0.001762778, 0.00521704, 0.003, 0.001];

% Place here the coefficient to adjust the points given through the tests in the sphere
corr_sensor_factor = [1.028112, 1.031493, 1.031425, 1.031247, 1.033897, 1.034449, 1.035083, 1.033594, 1.2384]; % OSRAM Coefficients
%corr_sensor_factor = [1.338262461, 1.044222867, 0.8604542965, 0.803963411, 0.6857015252, 0.788397613, 0.9612483281, 1] % Coefficients I took using sphere data

% Special XYZ Calibration Matrix for  Golden Device ("Created" Reference Device from gathered 
% production-device data) - Based on Full, Production-Wide Device-Characteristic Dataset (not a 
% per-device, so gives good, but not perfect results - Supports production without per-unit calibration) 
y_correction = [0.01396, 0.16748, 0.23538, 1.4275, 1.8867, 1.142, 0.46497, -0.027018, -0.24468, -0.019926];

% Conversion constant for PFD
ConstPP = (1e-3)/(6.02214076e23*2.998e8*6.62607015e-34);

figure(1)
h = plot(sensor_waveLength, sensor_values, '-o');
grid on;
xlim([400 1000]);
%ylim([0 1000]);
title("Light Sensors");
ylabel("Magnitude");
xlabel("Wavelength [nm]");

sensor_waveLength = [415, 445, 480, 515, 555, 590, 630, 680, 750, 900];

% Reads the data from the table and separates it into vectors
matrixValues = xlsread("General_Spec_Corr_Matrix.xlsx");

while(1)
    value = readline(arduino);
    F1 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F2 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F3 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F4 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F5 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F6 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F7 = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    F8 = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    Clear = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    NIR = str2double(value)/(GAIN256 * tint_ms);
   
    sensor_values = [F1, F2, F3, F4, F5, F6, F7, F8, Clear, NIR];
    
    sensor_values_corr = (sensor_values - corr_sensor_offset).*corr_sensor_factor;
    sensor_values_corr = sensor_values_corr ./ 2;  % Just a correction test value. It's not a confirmed value
    
    iluminance = y_correction * sensor_values_corr' * 683
    
    % Test conversion to PFD using sensor_values_corr as radiant power
    PFD = (sum(sensor_waveLength.*sensor_values_corr))*ConstPP
    
    % Test of PPFD just removing the clear and NIR sensors
    PPFD = (sum(sensor_waveLength(1,1:8).*sensor_values_corr(1,1:8)))*ConstPP
    
    set(h, 'XData', sensor_waveLength, 'YData', sensor_values_corr);
    pause(1);
    clc;
end

%% Plot sensor curves
clc;
close all;
clear all;

% Reads the data from the table of curves and separates it into vectors
sensorCurves = xlsread("Typical_AS7341_filters.xlsx");
WaveLenght = sensorCurves(3:end, 1);
F1curve = sensorCurves(3:end, 2);
F2curve = sensorCurves(3:end, 3);
F3curve = sensorCurves(3:end, 4);
F4curve = sensorCurves(3:end, 5);
F5curve = sensorCurves(3:end, 6);
F6curve = sensorCurves(3:end, 7);
F7curve = sensorCurves(3:end, 8);
F8curve = sensorCurves(3:end, 9);
Clearcurve = sensorCurves(3:end, 10);

figure(2);
plot(WaveLenght, F1curve);
title("Characteristic curve of the sensor filters");
hold on;
plot(WaveLenght, F2curve);
plot(WaveLenght, F3curve);
plot(WaveLenght, F4curve);
plot(WaveLenght, F5curve);
plot(WaveLenght, F6curve);
plot(WaveLenght, F7curve);
plot(WaveLenght, F8curve);
plot(WaveLenght, Clearcurve);
grid on;
xlabel("Wavelength [nm]");
ylabel("Magnitude");
legend(["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "Clear"]);

%% General Spectral Correction Matrix (GSCM) per Channel - Step size 1nm 
clc;
close all;
clear all;

% Reads the data from the table of curves and separates it into vectors
matrixValues = xlsread("General_Spec_Corr_Matrix.xlsx");
WaveLenght = matrixValues(2:end, 1);
matrix_GSCM = matrixValues(2:end, 2:end);
F1_GSCM = matrixValues(2:end, 2);
F2_GSCM = matrixValues(2:end, 3);
F3_GSCM = matrixValues(2:end, 4);
F4_GSCM = matrixValues(2:end, 5);
F5_GSCM = matrixValues(2:end, 6);
F6_GSCM = matrixValues(2:end, 7);
F7_GSCM = matrixValues(2:end, 8);
F8_GSCM = matrixValues(2:end, 9);
Clear_GSCM = matrixValues(2:end, 10);
NIR_GSCM = matrixValues(2:end, 11);

% figure(3);
% plot(WaveLenght, F1_GSCM);
% title("General Spectral Correction Matrix (GSCM) per Channel");
% hold on;
% plot(WaveLenght, F2_GSCM);
% plot(WaveLenght, F3_GSCM);
% plot(WaveLenght, F4_GSCM);
% plot(WaveLenght, F5_GSCM);
% plot(WaveLenght, F6_GSCM);
% plot(WaveLenght, F7_GSCM);
% plot(WaveLenght, F8_GSCM);
% plot(WaveLenght, Clear_GSCM);
% grid on;
% xlabel("Comprimento de onda [nm]");
% ylabel("Magnitude");

%% Spectral Reconstruction based on Channel data's
F1 = 0.011057;
F2 = 0.019664;
F3 = 0.028302;
F4 = 0.032492;
F5 = 0.034778;
F6 = 0.034016;
F7 = 0.039658;
F8 = 0.043168;
Clear = 0.127734;
NIR = 0.051806;

WaveLenght = 380:1:1000;

sensorValuesMatrix = [F1, F2, F3, F4, F5, F6, F7, F8, Clear, NIR]';

y_correction = [0.01396, 0.16748, 0.23538, 1.4275, 1.8867, 1.142, 0.46497, -0.027018, -0.24468, -0.019926];

iluminance = y_correction * sensorValuesMatrix * 683

reconstructedSpectro = matrix_GSCM * sensorValuesMatrix;

reconstructedSpectroNorm = reconstructedSpectro./max(reconstructedSpectro);

figure(4)
plot(WaveLenght, reconstructedSpectro, 'lineWidth', 1);
title("Reconstructed Spectrum");
xlabel("Wavelength [nm]");
ylabel("Sensitivity");
grid on;
xlim([380 1000]);

figure(5)
plot(WaveLenght, reconstructedSpectroNorm, 'lineWidth', 1);
title("Reconstructed normalized spectrum");
xlabel("Wavelength [nm]");
ylabel("Sensitivity");
grid on;
xlim([380 1000]);

%% Read sensors and plot spectrum from calculation using the General Spectral Correction Matrix (GSCM)
clc; 
clear all;
arduino  = serialport("COM13", 9600);

ATIME = 0;
ASTEP = 65534;
GAIN512 = 512; % Para F1 - F6
GAIN256 = 256; % Para F7, F8, Clear

tint_ms = (ATIME + 1) * (ASTEP + 1) * (2.78/1000);

WaveLenght = 380:1:1000;
reconstructedSpectro = zeros(1,(1000-380+1));
sensor_values = zeros(1,9);

% Finds the indices of the values closest to 400nm and 700nm
% in order to calculate the PPF
for i=1:length(WaveLenght)
    if(ge(WaveLenght(1, i),400))
        ind_inf = i;
        break
    end
end

for i=flip(1:length(WaveLenght))
    if(le(WaveLenght(1, i),700))
        ind_sup = i;
        break
    end
end

% Reads the data from the GSCM table
matrixValues = xlsread("General_Spec_Corr_Matrix.xlsx");
matrix_GSCM = matrixValues(2:end, 2:end);

% Reads the CIE1931 table
CIE1931 = xlsread("CIE1931_XYZ.xlsx");
CIE1931_Y = CIE1931(:,3);


% With an offset , the basic counts will be reduced, Offset can be DARK CURRENTS or other 
% constant vlaus those will affects the measurements permantly,
corr_sensor_offset = [0.00196979, 0.00724927, 0.00319381, 0.001314659, 0.001468153, 0.001858105, 0.001762778, 0.00521704, 0.003, 0.001];

% Place here the coefficient to adjust the points given through the tests in the sphere
corr_sensor_factor = [1.028112, 1.031493, 1.031425, 1.031247, 1.033897, 1.034449, 1.035083, 1.033594, 1.2384, 1.269416]; % OSRAM Coefficients
%corr_sensor_factor = [1.338262461, 1.044222867, 0.8604542965, 0.803963411, 0.6857015252, 0.788397613, 0.9612483281, 1] % Coefficients I took using sphere data


% Conversion constant for PFD
ConstPP = (1e-3)/(6.02214076e23*2.998e8*6.62607015e-34);

figure(1)
h = plot(WaveLenght, reconstructedSpectro, 'b');
grid on;
xlim([380 1000]);
title("Reconstructed spectrum of measured light");
ylabel("Sensitivity");
xlabel("Wavelength [nm]");

buffer_sensor_values = zeros(100, 10);
iteracao = 1;

% Loop to read the sensor data and calculate the spectrum
while(1)
    value = readline(arduino);
    F1 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F2 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F3 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F4 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F5 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F6 = str2double(value)/(GAIN512 * tint_ms);
    
    value = readline(arduino);
    F7 = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    F8 = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    Clear = str2double(value)/(GAIN256 * tint_ms);
    
    value = readline(arduino);
    NIR = str2double(value)/(GAIN256 * tint_ms);
   
    sensor_values = [F1, F2, F3, F4, F5, F6, F7, F8, Clear, NIR];
    
    sensor_values_corr = (sensor_values - corr_sensor_offset).*corr_sensor_factor;
    %sensor_values_corr = sensor_values;
    
    reconstructedSpectro = matrix_GSCM * sensor_values_corr';
   
    reconstructedSpectro = reconstructedSpectro ./ 2;  % Just a correction test value. It's not a confirmed value
    
    iluminance = reconstructedSpectro(1:401,1)' * CIE1931_Y * 683;
    
    PFD = (sum(WaveLenght'.*reconstructedSpectro))*ConstPP;
    PPFD = (sum(WaveLenght(1, ind_inf:ind_sup)'.*reconstructedSpectro(ind_inf:ind_sup,1)))*ConstPP;
    
    set(h, 'XData', WaveLenght, 'YData', reconstructedSpectro);
    pause(1);
    clc;
    if iteracao < 100
        buffer_sensor_values(iteracao, :) = sensor_values;
        display(iteracao);
        iteracao = iteracao + 1;
    end
end