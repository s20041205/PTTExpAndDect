%% Shu-Tyng Last modified on: Apr, 21, 2017
% Function of SBP detection
% Thesis: HOLTEK sensor module
% 
function [ppg_vanish, ppg_appear, sbp_found] = sbpDect_v2(pkPPGcuff, locpkPPGcuff, pkPPGsen, locpkPPGsen, prPPGsen, CuffDC, fs);
%
ppi_cuff = 0;
ppi_change = [];
ppg_vanish = [];
ppg_appear = [];
sbp_found = [];
for i = 5:length(locpkPPGcuff)
    if (i == 5)
        ppi_cuff = locpkPPGcuff(i) - locpkPPGcuff(i-1);
    else
        % Find the point where PPGcuff disappear
        ppi_cuff = [ppi_cuff; locpkPPGcuff(i) - locpkPPGcuff(i-1)];
        if (ppi_cuff(end) > ppi_cuff(end-1)*1.15) || (ppi_cuff(end) < ppi_cuff(end-1)*0.85)
            for ct_sen = length(locpkPPGsen):-1:1
            % diff (locpkPPGcuff - locpkPPGsen) to make sure that pkPPGcuff is a peak
                if (abs(locpkPPGcuff(i)-locpkPPGsen(ct_sen))<30)
                    ppi_change = i;
                    sbp_found = locpkPPGcuff(i);
                    fprintf('SBP found (%d): %0.2f\n', locpkPPGcuff(i), CuffDC(locpkPPGcuff(i)));
                end
            end
        end
    end
end
end % End of function
