%% Shu-Tyng Last modified on: Apr, 27, 2017
% Rawdata classification
% Thesis data read
clear; clc; close all; % close all
%% Parameter initialization
ct_fig = 1; % Initial figure counter
plot_en = 1; % Plot if enabled. (2:)
fs = 200; % 50 or 200
rec = [];
rec_data = [];
record = [];

% For SBP decision
rec_pk = [];
rec_pkData = [];
rec_sbp = [];
rec_pkSBP = [];
sbp = [];
locSBP = [];
loc_pk_slopeMax = [];
loc_sbp = [];
pt_sbp = [];
%% File selection
fprintf('Sample rate = %d\n', fs);
drt = '.\Holtek_Rawdata\Thesis_Module Data\';
[filename, drt] = uigetfile([drt, '*.txt'], 'Select a recorded data.');
if isequal(filename,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(drt, filename)]);
end
dataroute = fullfile(drt, filename);
fprintf('Data select: %s\n', filename);
%% Read & manage rawdata
% Holtek MCU with Sensor Module + Cuff Module
% Read ECG, PPGcuff, PPGsen, CuffAC, CuffDC, Battery
if (fs == 50) % For files of 50Hz
    [time, rECG, PPG1, PPG2, Battery, CuffDC, CuffAC] = textread(dataroute,'%s%d%d%d%d%d%d','headerlines',0);
elseif (fs == 200) % For files of 200Hz
    [date, time, rECG, PPG1, PPG2, Battery, CuffDC, CuffAC, record] = textread(dataroute,'%s%s%d%d%d%d%d%d%d','headerlines',0);
else
end
% samplerate = 200; % 200Hz
samples = 1:1:length(time);
% ttime = length(time)/fs; % Total seconds
 
% pvalue = rawPPG(:,1); % PPG value
rECG = rECG; %* 0.19 * 10^(-6); % 1bit = 0.19uV
PPG1 = PPG1; %* 0.8 * 10^(-3); % 1bit = 0.8mV
PPG2 = PPG2; %* 0.8 * 10^(-3); % 1bit = 0.8mV
CuffDC =  (CuffDC - 50.0989) / 0.7157; % CuffDC adjustment
 
% Set PPGsen & PPGcuff
PPGsen = PPG2; % PPG_sensor
PPGcuff = PPG1; % PPG_cuff
%% PPG peak detection
[pkPPGcuff, locpkPPGcuff, troughPPGcuff, loctroughPPGcuff, prPPGcuff] = PPGpkdec_WH(PPGcuff, fs);
locpkPPGcuff = locpkPPGcuff-35; % Delay shift
for p = 1:length(locpkPPGcuff)
   pkPPGcuff(p) = PPGcuff(locpkPPGcuff(p));
end
 
[pkPPGsen, locpkPPGsen, troughPPGsen, loctroughPPGsen, prPPGsen] = PPGpkdec_WH(PPGsen, fs);
locpkPPGsen = locpkPPGsen-35; % Delay shift
for p = 1:length(locpkPPGsen)
   pkPPGsen(p) = PPGsen(locpkPPGsen(p));
end
%% ECG peak detection
% Filtering
% fECG = filter(ecgLPF,rECG);
% ECG = filter(ecgHPF,fECG);
% P&T method
[pkECG, locpkECG, ECGdelay] = pan_tompkin(rECG, fs, 0);
locpkECG = locpkECG + 1; % Delay shift
for n = 1:length(locpkECG)
   pkECG(n) = rECG(locpkECG(n));
end
pkECG = pkECG';
locpkECG = locpkECG';
 
%% SBP decision
% [loc_pk_slopeMax, loc_sbp, pt_sbp] = sbpDect(pkPPGcuff, locpkPPGcuff, locpkPPGsen, prPPGsen, CuffDC, fs);
% [loc_pk_slopeMax, loc_sbp, pt_sbp] = sbpDect_v2(pkPPGcuff, locpkPPGcuff, locpkPPGsen, prPPGsen, CuffDC, fs);
ppi_cuff = [];
ppi_change = [];
% ppg_vanish = [];
ppg_appear = [];
pt_sbp = [];
ppi_temp = [];
ppi_ct = 0;
ppiCt_lim = 5;
% Find pkPPGcuff after CuffDCmax: deflation start
ct = 1;
while locpkPPGcuff(ct) < 1000
    ct = ct + 1;
end
ppi_temp = locpkPPGcuff(ct) - locpkPPGcuff(ct-1); % save refrence ppi
for ct_slope = ct: length(pkPPGcuff(:,1))
    if (pkPPGcuff(ct_slope) > pkPPGcuff(ct_slope-1)*1.2) && (pkPPGcuff(ct_slope) > 50) && (CuffDC(locpkPPGcuff(ct_slope))>=135)
       if isempty(loc_pk_slopeMax)
           loc_pk_slopeMax = ct_slope;
           fprintf('Deflation start (%d): %0.2f\n', locpkPPGcuff(ct_slope), CuffDC(locpkPPGcuff(ct_slope)));
       end
    end
end

for i = loc_pk_slopeMax:length(locpkPPGcuff)
    ppi_cuff = locpkPPGcuff(i) - locpkPPGcuff(i-1);
    if isempty(pt_sbp)
        if (ppi_cuff >= ppi_temp*0.8) && (ppi_cuff <= ppi_temp*1.2)
            ppi_ct = ppi_ct + 1;
            ppg_appear = [ppg_appear; i]; 
            if (ppi_ct > ppiCt_lim) && (ppi_ct-ppiCt_lim >= ppg_appear(end)-ppg_appear(end-ppiCt_lim))                    
                for ct = 2:length(ppg_appear)
                    if isempty(ppi_change) && (ppg_appear(ct)-ppg_appear(ct-1) == 1)
                        ppi_change = [ppi_change; ppg_appear(ct-1)];

                    pt_sbp = locpkPPGcuff(ppi_change(end)-1);
                    fprintf('SBP found (%d): %0.2f\n', pt_sbp, CuffDC(pt_sbp));
                    end
                end               
            end
        end
    end
end

%% PTT calculation
[ptt] = PTTcalc(locpkECG, locpkPPGsen);

%% Find CuffDC & CuffAC peaks
[pkCuffDC_max, locpkCuffDC_max] = CuffDCpkdec(CuffDC);
[pkCuffAC, locpkCuffAC, troughCuffAC, loctroughCuffAC, prCuffAC] = PPGpkdec_WH(CuffAC, fs);
locpkCuffAC = locpkCuffAC - 24; % Delay shift
for p = 1:length(locpkCuffAC)
   CuffAC(p) = CuffAC(locpkCuffAC(p));
end
 
%% Find record flag
for i = 1:length(record)
    if (record(i) == 1)
       rec = [rec; i];
       rec_data = [rec_data; rECG(i), -PPGcuff(i), -PPGsen(i), Battery(i), CuffDC(i), CuffAC(i)];
    end
end
%% Find loc_pk_slopeMax flag: for deflation start
for i = 1:length(loc_pk_slopeMax)
   rec_pk = [rec_pk; locpkPPGcuff(loc_pk_slopeMax(i))];
   rec_pkData = [rec_pkData; rECG(rec_pk(i)), -PPGcuff(rec_pk(i)), -PPGsen(rec_pk(i)), Battery(rec_pk(i)), CuffDC(rec_pk(i)), CuffAC(rec_pk(i))];
end
%% Find loc_sbp: for ppg appears
for i = 1:length(loc_sbp)
   rec_sbp = [rec_sbp; locpkPPGcuff(loc_sbp(i))];
   rec_pkSBP = [rec_pkSBP; rECG(rec_sbp(i)), -PPGcuff(rec_sbp(i)), -PPGsen(rec_sbp(i)), Battery(rec_sbp(i)), CuffDC(rec_sbp(i)), CuffAC(rec_sbp(i))];
end
%% Find pt_sbp: for sbp decision
for i = 1:length(pt_sbp)
   sbp = [sbp; pt_sbp(i)];
   locSBP = [locSBP; rECG(sbp(i)), -PPGcuff(sbp(i)), -PPGsen(sbp(i)), Battery(sbp(i)), CuffDC(sbp(i)), CuffAC(sbp(i))];
end
%% Plotting
if (plot_en == 1)
    %% Fig1. Plot ECG, PPGcuff, PPGsen
    x_L = 0;
    x_H = inf;
    y_L = -inf;
    y_H = inf;
   
    figure(ct_fig),
    subplot(3,1,1), plot(samples, rECG), hold on, grid on,
    plot(locpkECG, pkECG, 'ro'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,1)-50 rec_data(f,1)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,1)-50 rec_pkData(f,1)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,1)-50 rec_pkSBP(f,1)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,1)-50 locSBP(f,1)+50], 'g-.', 'LineWidth',2)
    end      
    xlabel ('Time (Samples)'), ylabel ('Voltage (uV)'),...
    title ('Rawdata of ECG'), axis ([x_L, x_H, y_L, y_H]),
 
    subplot(3,1,2), plot(samples, -PPGcuff), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,2)-50 rec_data(f,2)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,2)-50 rec_pkData(f,2)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,2)-50 rec_pkSBP(f,2)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,2)-50 locSBP(f,2)+50], 'g-.', 'LineWidth',2)
    end      
    plot(locpkPPGcuff(:,1), -pkPPGcuff(:,1), 'ro'), plot(loctroughPPGcuff, -troughPPGcuff, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGcuff'), axis ([x_L, x_H, y_L, y_H]),
 
    subplot(3,1,3), plot(samples, -PPGsen), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,3)-50 rec_data(f,3)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,3)-50 rec_pkData(f,3)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,3)-50 rec_pkSBP(f,3)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,3)-50 locSBP(f,3)+50], 'g-.', 'LineWidth',2)
    end
    plot(locpkPPGsen(:,1), -pkPPGsen(:,1), 'ro'), plot(loctroughPPGsen, -troughPPGsen, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGsensor'), axis ([x_L, x_H, y_L, y_H]),
    ct_fig = ct_fig +1;
  
% elseif (plot_en == 2)
    %% Fig2. Plot CuffDC, PPGcuff, PPGsen
    x_L = 1000;
    x_H = inf;
    y_L = 0;
    y_H = inf;  
    figure(ct_fig),
    subplot(4,1,1), plot(samples, CuffDC), hold on, grid on,
    plot(locpkCuffDC_max, pkCuffDC_max, 'gx'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,5)-50 rec_data(f,5)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,5)-50 rec_pkData(f,5)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,5)-50 rec_pkSBP(f,5)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,5)-50 locSBP(f,5)+50], 'g-.', 'LineWidth',2)
    end
    xlabel ('Time (Samples)'), ylabel ('Pressure (mmHg)'),...
    title ('CuffDC'), axis ([x_L, x_H, y_L, y_H]),
 
    y_L = 1200;
    y_H = 2800;
    subplot(4,1,2), plot(samples, CuffAC), hold on, grid on,
    plot(locpkCuffAC, pkCuffAC, 'ro'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,6)-50 rec_data(f,6)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,6)-50 rec_pkData(f,6)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,6)-50 rec_pkSBP(f,6)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,6)-50 locSBP(f,6)+50], 'g-.', 'LineWidth',2)
    end
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('CuffAC'), axis ([x_L, x_H, y_L, y_H]),
 
    y_L = -inf;
    y_H = inf;
    subplot(4,1,3), plot(samples, -PPGcuff), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,2)-50 rec_data(f,2)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,2)-50 rec_pkData(f,2)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,2)-50 rec_pkSBP(f,2)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,2)-50 locSBP(f,2)+50], 'g-.', 'LineWidth',2)
    end
    plot(locpkPPGcuff(:,1), -pkPPGcuff(:,1), 'ro'), plot(loctroughPPGcuff, -troughPPGcuff, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGcuff'), axis ([x_L, x_H, y_L, y_H]),
 
    y_L = -inf;
    y_H = inf;
    subplot(4,1,4), plot(samples, -PPGsen), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,3)-50 rec_data(f,3)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,3)-50 rec_pkData(f,3)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,3)-50 rec_pkSBP(f,3)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,3)-50 locSBP(f,3)+50], 'g-.', 'LineWidth',2)
    end
    plot(locpkPPGsen(:,1), -pkPPGsen(:,1), 'ro'), plot(loctroughPPGcuff, -troughPPGsen, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGsensor'), axis ([x_L, x_H, y_L, y_H]),
    ct_fig = ct_fig +1;
  

    %% Fig3. Plot CuffDC, CuffAC, Battery 
    figure(ct_fig),
    subplot(3,1,1), hold on, grid on,
    [hAx,hLine1,hLine2] = plotyy(samples, CuffDC, samples, CuffAC);
    plot(locpkCuffDC_max, pkCuffDC_max, 'ro'), plot(locpkCuffAC, pkCuffAC, 'ro'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,5)-50 rec_data(f,5)+50], 'm')
    end
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,6)-50 rec_data(f,6)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,5)-50 rec_pkData(f,5)+50], 'm--')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,6)-50 rec_pkData(f,6)+50], 'm--')
    end      
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,5)-50 rec_pkSBP(f,5)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,6)-50 rec_pkSBP(f,6)+50], 'm:', 'LineWidth',2)
    end
   for f = 1:length(sbp)
       plot([sbp(f,:) sbp(f,:)], [locSBP(f,5)-50 locSBP(f,5)+50], 'g-.', 'LineWidth',2)
   end
   for f = 1:length(sbp)
       plot([sbp(f,:) sbp(f,:)], [locSBP(f,6)-50 locSBP(f,6)+50], 'g-.', 'LineWidth',2)
   end
   xlabel ('Time (Samples = 200Hz)'),
   ylabel (hAx(1),'Cuff pressure (mmHg)'), ylabel (hAx(2),'ADC value'),
   title ('CuffDC & CuffAC'),
  
   x_L = -inf;
   x_H = inf;
   y_L = -inf;
   y_H = inf;
   subplot(3,1,2), plot(samples, Battery), hold on, grid on,
   for f = 1:length(rec)
       plot([rec(f,:) rec(f,:)], [rec_data(f,4)-50 rec_data(f,4)+50], 'm')
   end
   for f = 1:length(rec_pk)
       plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,4)-50 rec_pkData(f,4)+50], 'm--')
   end
   for f = 1:length(rec_sbp)
       plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,4)-50 rec_pkSBP(f,4)+50], 'm:', 'LineWidth',2)
   end
   for f = 1:length(sbp)
       plot([sbp(f,:) sbp(f,:)], [locSBP(f,4)-50 locSBP(f,4)+50], 'g-.', 'LineWidth',2)
   end
   xlabel ('Time (Samples)'), ylabel ('ADC value'),...
       title ('Battery'), axis ([x_L, x_H, y_L, y_H]),
   
   subplot(313)
   
   ct_fig = ct_fig +1;
elseif (plot_en == 2)
       %% Fig4. Plot CuffDC & CuffAC, ECG, PPGcuff, PPGsensor 
    figure(ct_fig),
    subplot(4,1,1), hold on, grid on,
    [hAx,hLine1,hLine2] = plotyy(samples, CuffDC, samples, CuffAC);
    plot(locpkCuffDC_max, pkCuffDC_max, 'ro'), plot(locpkCuffAC, pkCuffAC, 'ro'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,5)-50 rec_data(f,5)+50], 'm')
    end
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,6)-50 rec_data(f,6)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,5)-50 rec_pkData(f,5)+50], 'm--')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,6)-50 rec_pkData(f,6)+50], 'm--')
    end      
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,5)-50 rec_pkSBP(f,5)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,6)-50 rec_pkSBP(f,6)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
       plot([sbp(f,:) sbp(f,:)], [locSBP(f,5)-50 locSBP(f,5)+50], 'g-.', 'LineWidth',2)
    end
    for f = 1:length(sbp)
       plot([sbp(f,:) sbp(f,:)], [locSBP(f,6)-50 locSBP(f,6)+50], 'g-.', 'LineWidth',2)
    end
    xlabel ('Time (Samples = 200Hz)'),
    ylabel (hAx(1),'Cuff pressure (mmHg)'), ylabel (hAx(2),'ADC value'),
    title ('CuffDC & CuffAC'),
    
    x_L = -inf;
    x_H = inf;
    y_L = -inf;
    y_H = inf;
    subplot(4,1,2), plot(samples, rECG), hold on, grid on,
    plot(locpkECG, pkECG, 'ro'),
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,1)-50 rec_data(f,1)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,1)-50 rec_pkData(f,1)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,1)-50 rec_pkSBP(f,1)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,1)-50 locSBP(f,1)+50], 'g-.', 'LineWidth',2)
    end      
    xlabel ('Time (Samples)'), ylabel ('Voltage (uV)'),...
    title ('Rawdata of ECG'), axis ([x_L, x_H, y_L, y_H]),
 
    subplot(4,1,3), plot(samples, -PPGcuff), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,2)-50 rec_data(f,2)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,2)-50 rec_pkData(f,2)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,2)-50 rec_pkSBP(f,2)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,2)-50 locSBP(f,2)+50], 'g-.', 'LineWidth',2)
    end      
    plot(locpkPPGcuff(:,1), -pkPPGcuff(:,1), 'ro'), plot(loctroughPPGcuff, -troughPPGcuff, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGcuff'), axis ([x_L, x_H, y_L, y_H]),
 
    subplot(4,1,4), plot(samples, -PPGsen), hold on, grid on,
    for f = 1:length(rec)
        plot([rec(f,:) rec(f,:)], [rec_data(f,3)-50 rec_data(f,3)+50], 'm')
    end
    for f = 1:length(rec_pk)
        plot([rec_pk(f,:) rec_pk(f,:)], [rec_pkData(f,3)-50 rec_pkData(f,3)+50], 'm--')
    end
    for f = 1:length(rec_sbp)
        plot([rec_sbp(f,:) rec_sbp(f,:)], [rec_pkSBP(f,3)-50 rec_pkSBP(f,3)+50], 'm:', 'LineWidth',2)
    end
    for f = 1:length(sbp)
        plot([sbp(f,:) sbp(f,:)], [locSBP(f,3)-50 locSBP(f,3)+50], 'g-.', 'LineWidth',2)
    end
    plot(locpkPPGsen(:,1), -pkPPGsen(:,1), 'ro'), plot(loctroughPPGsen, -troughPPGsen, 'gx'),
    xlabel ('Time (Samples)'), ylabel ('ADC value'),...
    title ('Rawdata of PPGsensor'), axis ([x_L, x_H, y_L, y_H]),
    ct_fig = ct_fig +1;
 
else
end
