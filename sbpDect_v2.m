%% Shu-Tyng Last modified on: Apr, 21, 2017
% Function of SBP detection
% Thesis: HOLTEK sensor module
% 
function [ppg_vanish, ppg_appear, sbp_found] = sbpDect_v2(pkPPGcuff, locpkPPGcuff, pkPPGsen, locpkPPGsen, prPPGsen, CuffDC, fs);
%
% ppi_cuff = 0;
% ppi_change = [];
% ppg_vanish = [];
% ppg_appear = [];
% sbp_found = [];
% for i = 5:length(locpkPPGcuff)
%     if (i == 5)
%         ppi_cuff = locpkPPGcuff(i) - locpkPPGcuff(i-1);
%     else
%         % Find the point where PPGcuff disappear
%         ppi_cuff = [ppi_cuff; locpkPPGcuff(i) - locpkPPGcuff(i-1)];
%         if (ppi_cuff(end) > ppi_cuff(end-1)*1.15) || (ppi_cuff(end) < ppi_cuff(end-1)*0.85)
%             for ct_sen = length(locpkPPGsen):-1:1
%             % diff (locpkPPGcuff - locpkPPGsen) to make sure that pkPPGcuff is a peak
%                 if (abs(locpkPPGcuff(i)-locpkPPGsen(ct_sen))<30)
%                     ppi_change = i;
%                     sbp_found = locpkPPGcuff(i);
%                     fprintf('SBP found (%d): %0.2f\n', locpkPPGcuff(i), CuffDC(locpkPPGcuff(i)));
%                 end
%             end
%         end
%     end
% end
% end % End of function
ppi_cuff = [];
ppi_change = [];
% ppg_vanish = [];
ppg_appear = [];
pt_sbp = [];
ppi_temp = [];
ppi_ct = 0;
ppiCt_lim = 10;
% Find pkPPGcuff after CuffDCmax: deflation start
ct = 1;
while locpkPPGcuff(ct) < 1000
    ct = ct + 1;
end
ppi_temp = locpkPPGcuff(ct) - locpkPPGcuff(ct-1); % save refrence ppi
% Start point
for ct_slope = ct: length(pkPPGcuff(:,1))
    if ((locpkPPGcuff(ct_slope)-locpkPPGcuff(ct_slope-1)) >= ppi_temp*1.5) && (pkPPGcuff(ct_slope) > 50) %&& (CuffDC(locpkPPGcuff(ct_slope))>=135)
       if isempty(loc_pk_slopeMax)
           loc_pk_slopeMax = ct_slope;
           fprintf('Deflation start (%d): %0.2f\n', locpkPPGcuff(ct_slope), CuffDC(locpkPPGcuff(ct_slope)));
       end
    end
end

for i = loc_pk_slopeMax:length(locpkPPGcuff)
    ppi_cuff = locpkPPGcuff(i) - locpkPPGcuff(i-1);
    if isempty(pt_sbp)
        if (ppi_cuff >= ppi_temp*0.8) && (ppi_cuff <= ppi_temp*1.2) && (CuffDC(locpkPPGcuff(i))<= CuffDC(locpkPPGcuff(loc_pk_slopeMax)))
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