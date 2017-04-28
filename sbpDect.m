%% Shu-Tyng Last modified on: Apr, 21, 2017
% Function of SBP detection
% Thesis: HOLTEK sensor module
% 
function [loc_pk_slopeMax, loc_sbp, pt_sbp] = sbpDect(pkPPGcuff, locpkPPGcuff, locpkPPGsen, prPPGsen, CuffDC, fs);
%% Initialization
loc_pk_slopeMax = [];
loc_sbp = [];
pt_sbp = [];
%% i. ppi check
ppi_ThrL = 0.7;
ppi_ThrH = 1.3;
% PPGcuff = PPG_cuff
locpkPPGcuff(:,2) = 0; % % Initialized for peak_slope storage
for ppi_ck = 2:length(locpkPPGcuff(:,1))
    locpkPPGcuff(ppi_ck, 2) = locpkPPGcuff(ppi_ck) - locpkPPGcuff(ppi_ck-1);
    if (locpkPPGcuff(ppi_ck, 2)<=40) || (locpkPPGcuff(ppi_ck, 2)>=400)
        locpkPPGcuff(ppi_ck, 2) = -1; % -1 for unreasonable
    end  
    % PPI[n-1]*0.7 < PPI[n] < PPI[n-1]*1.3
    if ((locpkPPGcuff(ppi_ck, 2)<=locpkPPGcuff(ppi_ck-1, 2)*0.7)&&(locpkPPGcuff(ppi_ck-1, 2)~=-1)) || ((locpkPPGcuff(ppi_ck, 2)>=locpkPPGcuff(ppi_ck-1, 2)*1.3)&&(locpkPPGcuff(ppi_ck-1, 2)~=-1))
        if (locpkPPGcuff(ppi_ck-1, 2)==-1)
            locpkPPGcuff(ppi_ck, 2) = locpkPPGcuff(ppi_ck, 2);
        else
            locpkPPGcuff(ppi_ck, 2) = -1; % -1 for unreasonable
        end  
    end
end
% PPGsen = PPG_sensor
locpkPPGsen(:,2) = 0; % Initialized for peak_slope storage
for ppi_ck = 2:length(locpkPPGsen(:,1))
    locpkPPGsen(ppi_ck, 2) = locpkPPGsen(ppi_ck) - locpkPPGsen(ppi_ck-1);
    if (locpkPPGsen(ppi_ck, 2)<=40) || (locpkPPGsen(ppi_ck, 2)>=400)
        locpkPPGsen(ppi_ck, 2) = -1; % -1 for unreasonable
    end  
    % PPI[n-1]*0.7 < PPI[n] < PPI[n-1]*1.3
    if ((locpkPPGsen(ppi_ck, 2)<=locpkPPGsen(ppi_ck-1, 2)*0.7)&&(locpkPPGsen(ppi_ck-1, 2)~=-1)) || ((locpkPPGsen(ppi_ck, 2)>=locpkPPGsen(ppi_ck-1, 2)*1.3)&&(locpkPPGsen(ppi_ck-1, 2)~=-1))
        if (locpkPPGsen(ppi_ck-1, 2)==-1)
            locpkPPGsen(ppi_ck, 2) = locpkPPGsen(ppi_ck, 2);
        else
            locpkPPGsen(ppi_ck, 2) = -1; % -1 for unreasonable
        end  
    end
end
%% ii. pk slope check
pkPPGcuff(:,2) = 0; % Initialized for slope storage
for ct_slope = 2:length(pkPPGcuff(:,1))
    pkPPGcuff(ct_slope,2) = abs(pkPPGcuff(ct_slope, 1) - pkPPGcuff(ct_slope-1, 1));
end
%% iii. SBP decision
% Step1. ppg dissepar => Point of cuff inflation stopped and deflation started

% Finding sudden changes in pkPPGcuff(:,2)
% pk_slopeMax = [pkPPGcuff, slope_pkPPGcuff]
ct = 1;
while locpkPPGcuff(ct) < 1000
    ct = ct + 1;
end
for ct_slope = ct: length(pkPPGcuff(:,1))
    if (pkPPGcuff(ct_slope,1) > pkPPGcuff(ct_slope-1, 1)*1.2) && (pkPPGcuff(ct_slope,1) > 50) && (CuffDC(locpkPPGcuff(ct_slope))>=135)
       loc_pk_slopeMax = [loc_pk_slopeMax; ct_slope];
       fprintf('Deflation start (%d): %0.2f\n', locpkPPGcuff(ct_slope), CuffDC(locpkPPGcuff(ct_slope)));
    end
end
 
% Step2. Find the point when ppg appears
% Find the point where slope begins to rise suddenly

for ct_slope = loc_pk_slopeMax(end): length(pkPPGcuff(:,1))
    if isempty(loc_sbp)
        if (pkPPGcuff(ct_slope,1) < pkPPGcuff(ct_slope-1,1)*0.8)
            loc_sbp = ct_slope;
            fprintf('Start of SBP finding = CuffDC(%d): %0.2f\n', locpkPPGcuff(ct_slope), CuffDC(locpkPPGcuff(ct_slope)));
        end
    end
end
 
% Step3. Find SBP
% Search-back where the point begins periodically

if isempty(loc_pk_slopeMax)
    loc_pk_slopeMax = 5;
end
if isempty(loc_sbp)
    loc_sbp = length(pkPPGcuff)-5;
end
t = loc_pk_slopeMax;
def_period = CuffDC(locpkPPGcuff(loc_pk_slopeMax))-CuffDC(locpkPPGcuff(loc_sbp))
if (def_period <= 70)
    def_lim = 1000
elseif (def_period > 70)
    def_lim = 2000
end

while (locpkPPGcuff(t) < locpkPPGcuff(loc_pk_slopeMax)+def_lim)
    t = t + 1;
end
fprintf('End of SBP finding = CuffDC(%d): %0.2f\n', locpkPPGcuff(t), CuffDC(locpkPPGcuff(t)));
avgPPI = 60*fs/mean(prPPGsen);
fprintf('Average PR: %0.2f\n', mean(prPPGsen));
for ct = loc_sbp:-1:t % PPGcuff
    for ct_sen = length(locpkPPGsen):-1:t % PPGsen
%         if isempty(pt_sbp)
            if (abs(locpkPPGcuff(ct,1)-locpkPPGsen(ct_sen,1))<30) && ((pkPPGcuff(ct,2)-pkPPGcuff(ct-1,2)>5) || pkPPGcuff(ct+1,2)-pkPPGcuff(ct,2)>1) && (avgPPI*ppi_ThrL<=locpkPPGsen(ct_sen,2)<=avgPPI*ppi_ThrH)
                % 20(best), 25
                pt_sbp = locpkPPGsen(ct_sen-1);
                fprintf('SBP found (%d): %0.2f\n', locpkPPGsen(ct_sen-1), CuffDC(locpkPPGsen(ct_sen-1)));
%             end
        end
    end
end
end
