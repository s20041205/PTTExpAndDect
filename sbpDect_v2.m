%% Shu-Tyng Last modified on: May, 05, 2017
% Function of SBP detection
% Thesis: HOLTEK sensor module
% 
function [loc_pk_slopeMax, loc_sbp, pt_sbp] = sbpDect_v2(pkPPGcuff, locpkPPGcuff, locpkPPGsen, CuffDC, fs);
ppi_cuff = [];
ppi_stable = [];
ppg_appear = [];
pt_sbp = [];
% pr_sbp = [];
ppi_temp = [];
ppi_ct = 0;
ppiCt_lim = 6;
ct_sen = 1;
loc_sbp = [];
loc_pk_slopeMax = [];
% Find pkPPGcuff after CuffDCmax: deflation start
ct = 1;
while locpkPPGcuff(ct) < 1000
    ct = ct + 1;
end
ppi_temp = locpkPPGcuff(ct) - locpkPPGcuff(ct-1); % save refrence ppi
pr_temp = 60*fs/ppi_temp;
fprintf('PR: %0.2f bpm (w/ ppi = %d)\n', pr_temp, ppi_temp);
% Start point: deflation
for ct_slope = ct: length(pkPPGcuff(:,1))
    if (((locpkPPGcuff(ct_slope)-locpkPPGcuff(ct_slope-1)) >= ppi_temp*1.2) || (locpkPPGcuff(ct_slope)-locpkPPGcuff(ct_slope-1)) <= ppi_temp*0.8) && (pkPPGcuff(ct_slope) > 50) %&& (CuffDC(locpkPPGcuff(ct_slope))>=135)
       if isempty(loc_pk_slopeMax)
           loc_pk_slopeMax = ct_slope;
           fprintf('Deflation start (%d): %0.2f\n', locpkPPGcuff(ct_slope), CuffDC(locpkPPGcuff(ct_slope)));
       end
    end
end
% Periodic
for i = loc_pk_slopeMax:length(locpkPPGcuff)
    ppi_cuff = locpkPPGcuff(i) - locpkPPGcuff(i-1);
    if isempty(pt_sbp)
        if (ppi_cuff>=ppi_temp*0.8) && (ppi_cuff<=ppi_temp*1.2) && (CuffDC(locpkPPGcuff(i))<= CuffDC(locpkPPGcuff(loc_pk_slopeMax)))
            ppi_ct = ppi_ct + 1;
            ppg_appear = [ppg_appear; i];
            if (length(ppg_appear)>=2) && (ppg_appear(end)-ppg_appear(end-1)<=2)
                ppi_stable = [ppi_stable; ppg_appear(end-1)];
                if (length(ppi_stable)>=ppiCt_lim) && (ppi_stable(ppiCt_lim)-ppi_stable(1) <= ppiCt_lim+4)                    
                    for ct_sen = 1:length(locpkPPGsen)
                        if(abs(locpkPPGcuff(ppi_stable(1)-1)-locpkPPGsen(ct_sen))< 50)
                            pt_sbp = locpkPPGsen(ct_sen);
                            fprintf('SBP found (%d): %0.2f\n', pt_sbp, CuffDC(pt_sbp));                         
                        end
                    end
                elseif (length(ppi_stable)>=ppiCt_lim) && (ppi_stable(ppiCt_lim)-ppi_stable(1) > ppiCt_lim+4)
                    % move ppiCt_lim to the next point fit the requested range
                    for ct_mw = 1:length(ppi_stable)-ppiCt_lim
                        if (ppi_stable(ct_mw+ppiCt_lim)-ppi_stable(ct_mw) <= ppiCt_lim+4)
                            for ct_sen = 1:length(locpkPPGsen)
                                if(abs(locpkPPGcuff(ppi_stable(ct_mw)-1)-locpkPPGsen(ct_sen))< 50)
                                    pt_sbp = locpkPPGsen(ct_sen);
                                    fprintf('SBP found (%d): %0.2f\n', pt_sbp, CuffDC(pt_sbp));
                                end
                            end                            
                        end
                    end
                end
            end
        end
    end
end % End of function
