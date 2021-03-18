%% Carrier Smoothing:



prob_smooth = find_carrier_phase(prob_sort);

prob_smooth(1).rho_smooth = [];
prob_smooth = carrier_smoothing(prob_smooth);
%%
prob_sort(1).rho_smooth = 0;

for ind_t = 1:length(prob_smooth)
    for ind_sat = 1:length(prob_smooth(ind_t).Svid)
        prob_sort(ind_t).rho_meas(ind_sat) = prob_smooth(ind_t).rho_smooth(ind_sat);
    end
end

%% Plot LLA:
x_lla = zeros(length(prob_smooth),3); cmat = parula(length(prob_smooth));
for ind=1:length(prob_smooth)
    x_lla(ind,:) = ecef2lla(prob_smooth(ind).x0_kf','WGS84'); 
end
figure;
for i=1:length(prob_smooth)
    geoplot(x_lla(i,1), x_lla(i,2),'.','Color',cmat(i,:));hold on;
end
title('Receiver position in LLA frame using WGS84 datum');
colorbar('southoutside','TickLabelInterpreter','latex','FontSize',24,...
    'TicksMode','manual','Ticks',[0, 1], 'TickLabels',{'$t = 0$', '$t = t_{end}$'})

%% Functions:
function prob_smooth = find_carrier_phase(prob_sort)
    prob_smooth = struct();
    c = 299792458;
    for ind_t = 1:length(prob_sort)
        for ind_sat = 1:length(prob_sort(ind_t).Svid)
            prob_smooth(ind_t).utcTimeMillis(ind_sat) = prob_sort(ind_t).utcTimeMillis(ind_sat);
            prob_smooth(ind_t).Svid(ind_sat) = prob_sort(ind_t).Svid(ind_sat);
            prob_smooth(ind_t).rho_meas(ind_sat) = prob_sort(ind_t).rho_meas(ind_sat);
            prob_smooth(ind_t).Adrm(ind_sat) = prob_sort(ind_t).AccumulatedDeltaRangeMeters(ind_sat);
            lam = c/prob_sort(ind_t).CarrierFrequencyHz(ind_sat);
            prob_smooth(ind_t).phi(ind_sat) = prob_sort(ind_t).AccumulatedDeltaRangeMeters(ind_sat)/-lam;                
        end
    end
end

% WRONG!!!!!
function prob_smooth = carrier_smoothing(prob_smooth)
    len_w = 20; % Measurements
    for ind_t = 1:length(prob_smooth)
        for ind_sat = 1:length(prob_smooth(ind_t).Svid)
            if ind_t == 1
                rho_prevs = [prob_smooth(1:ind_t).rho_meas];
                del_phi_t_i = 0;
            else
                rho_prevs = [prob_smooth(1:ind_t-1).rho_meas];
                del_phi_t_i = prob_smooth(ind_t).Adrm(ind_sat);
            end
            
            M = length(rho_prevs);
            M_opt = min(M, len_w);
            rho_ti = prob_smooth(ind_t).rho_meas(ind_sat);
            rho_bar_tm = mean(rho_prevs(end:-1:M_opt));
            rho_bar_ti = 1/M_opt*rho_ti + (M_opt-1)/M_opt*(rho_bar_tm + del_phi_t_i);
            prob_smooth(ind_t).rho_smooth(ind_sat) = rho_bar_ti;
        end
    end
end

