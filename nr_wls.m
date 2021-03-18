%{
Filename: basic_solution.m
by, Shivam Soni - 03/12/2021
1) Reads the gnss log file "gnss_log.csv",
2) Reads the hourly gps ephemeris file "ephem_file.csv",
3) Performs basic localization algorithm: Newton-Raphson
%}

clearvars; clc; close all; 
set(0,'defaultTextInterpreter','latex');set(0,'defaultAxesFontSize',20);
set(0,'defaultLegendInterpreter','latex');

%% Initialize and find satellite positions:

% Extract the relavent Ephemeris data from the given ephemeris file:
ephem_file = 'ephem_file.csv'; % Ephemeris filename
prob_file = 'gnss_log.csv'; % gnss_log filename

% Create problem structure:
prob = create_prob(prob_file);

% Process the hourly ephemeris data
ephem_struct = ephem_process(ephem_file,prob);

%% SV clock Biases:
ephem_struct(1).B = [];
B = GpsEph2Dtsv(ephem_struct, [ephem_struct.Toe]);

for eph_ind = 1:length(ephem_struct)
    ephem_struct(eph_ind).B = B(eph_ind);
end

%% Find Satellite positions:
len_gnss = length(prob);
for ind = 1:len_gnss
    prob(ind) = ephem2pos(ephem_struct, prob(ind));
end

% Sort the problem according to time:
prob_sort = time_sort(prob);
%% Iono-Free corrections:
prob_sort(1).rho_meas_IF = [];
prob_sort(1).Svid_unique = [];
% Compute Iono-Free pseudoranges:
T_tot = length(prob_sort);
for ind_t = 1:T_tot
    prob_sort(ind_t) = rho_iono_free(prob_sort(ind_t),'meas');
end

prob_sort(1).sat_pos_unique = [];
% Reduce the Satellite Positions for iono free:
for ind_t = 1:T_tot
    prob_sort(ind_t) = sat_pos_reduction(prob_sort(ind_t));
end

%% Acceleration, Gyro, and Mag:
accel = table2array(readtable('accel.csv','ReadVariableNames',0));
gyro = table2array(readtable('gyro.csv','ReadVariableNames',0));
mag = table2array(readtable('mag.csv','ReadVariableNames',0));

prob_sort(1).a = []; % Accel
prob_sort(1).w = []; % Gyro
prob_sort(1).mag = []; % Magnetometer

prob_sort = gyro_accel_mag_sort(accel, gyro, mag, prob_sort);

%% Calibrate Accel, Gyro, and Mag:
prob_sort(1).a_cal = [];
prob_sort = accel_calib(prob_sort);

prob_sort(1).w_cal = [];
prob_sort = gyro_calib(prob_sort);

prob_sort(1).mag_cal = [];
prob_sort = mag_calib(prob_sort);

%% Newton-Raphson:
% Initial Guess
prob_sort(1).x0 = [0;0;0]; % Receiver Position
prob_sort(1).bu = 0; % Receiver Clock Bias
prob_sort(1).rho_theo = [];
prob_sort(1).rho_theo_IF = [];

% Update the position estiamate over each timestep:
T_tot = length(prob_sort);
for ind_t = 1:T_tot % Since last time instant is with no observations
    if length(unique(prob_sort(ind_t).Svid)) < 4 % Ignore if less than 4 unique SV's
        % Save the previous position an time estimates for the next timestep
        prob_sort(ind_t+1).x0 = prob_sort(ind_t-1).x0; 
        prob_sort(ind_t+1).bu = prob_sort(ind_t-1).bu;
        continue
    end
    prob_sort(ind_t) = NR_GPS(prob_sort(ind_t)); % Perform Newton-Raphson Algorithm
    % Use current position to guide the next timestep's initial guess
    if ind_t ~= T_tot % Until the end of the time range
        prob_sort(ind_t+1).x0 = prob_sort(ind_t).x0;
        prob_sort(ind_t+1).bu = prob_sort(ind_t).bu;
    end
end
%%
% Sanity check:
x_lla = zeros(length(prob_sort),3); cmat = parula(length(prob_sort));
for ind=1:T_tot
    X0 = prob_sort(1).x0'; 
    Xi = prob_sort(ind).x0'; 
    if norm(Xi - X0) < 50e3 % 50 km
        x_lla(ind,:) = ecef2lla(prob_sort(ind).x0','WGS84'); 
    end
end
SZ = size(x_lla);
% Remove zeros:
x_lla = reshape(x_lla, SZ(1)*SZ(2), 1);
x_lla(x_lla == 0) = [];
x_lla = reshape(x_lla, [], 3);
%%
figure;
for i=1:size(x_lla, 1)
    geoplot(x_lla(i,1), x_lla(i,2),'.','Color',cmat(i,:));hold on;
end
title('Receiver position in LLA frame using WGS84 datum');
colorbar('southoutside','TickLabelInterpreter','latex','FontSize',24,...
    'TicksMode','manual','Ticks',[0, 1], 'TickLabels',{'$t = 0$', '$t = t_{end}$'})

%% Save LLA:

writematrix(x_lla, 'LLA.csv')
%% Functions:

function ephem_prob = ephem_process(ephem_file, prob)
    ephem_struct = table2struct(readtable(ephem_file,'ReadVariableNames',1));

    t_prob_ind = 1;
    while prob(1).ReceivedSvTimeNanos/1e9 - ephem_struct(t_prob_ind).Toe > 0
        t_prob_ind = t_prob_ind+1;
    end
    
    TOE_prob = ephem_struct(t_prob_ind-1).Toe;
    ind = t_prob_ind-1;
    ephem_prob = ephem_struct(1); % Create a struct and copy all the fields:

    while ephem_struct(ind).Toe == TOE_prob
        ephem_prob(t_prob_ind - ind) = ephem_struct(ind);
        ind = ind-1;
    end
    e_ind = ind;
    ephem_prob = flip(ephem_prob');
    
    for ind = 1:32
        if ephem_prob(ind).PRN ~= ind
            ephem_prob = [ephem_prob(1:ind-1); find_eph(ind, ephem_struct, e_ind);
                ephem_prob(ind:end)];
        end
    end
    
    ind = 1;
    len_ephem_prob = length(ephem_prob);
    while ind <= len_ephem_prob
         if ephem_prob(ind).PRN == 100
             ephem_prob(ind) = [];
             len_ephem_prob = length(ephem_prob);
         end
         ind = ind + 1;
    end
    function eph_out = find_eph(PRN, ephem_struct, e_ind)
        temp_ind = e_ind;
        while temp_ind > 0 && ephem_struct(temp_ind).PRN ~= PRN 
            temp_ind = temp_ind - 1;
        end
        if temp_ind <= 0 
            eph_out = ephem_struct(1);
            eph_out.PRN = 100;
        else
            eph_out = ephem_struct(temp_ind);
        end
    end
end

function prob = create_prob(prob_file)
    prob = table2struct(readtable(prob_file,'ReadVariableNames',1));
    prob(1).sat_pos_calc = []; % ECEF satellite Position Vector
    prob(1).t_k = []; % t_k
    prob(1).B = []; % B
end

function prob = ephem2pos(ephem_struct, prob)
    c = 299792458; prn_ind = find([ephem_struct.PRN] == prob.Svid);
    if ~isempty(prn_ind)
        prob.t_k = prob.ReceivedSvTimeNanos/1e9 - ephem_struct(prn_ind).Toe;
        [~, dtsvS, ~, dtsvSdot] = GpsEph2Pvt(ephem_struct(prn_ind), ...
            [ephem_struct(prn_ind).GPS_Week, prob.ReceivedSvTimeNanos/1e9]);
        prob.B = dtsvS*c + dtsvSdot*c*prob.t_k;
        prob.sat_pos_calc = ephem2cart(ephem_struct(prn_ind), prob);
    else
        error('PRN %d not found in the given ephemeris data', prob.PRNS);
    end
end

function sat_pos_calc = ephem2cart(ephem_struct, prob)
    mu = 3.986005e14; O_dot_E = 7.2921151467e-5;
    a = ephem_struct.Asqrt^2;
    n = sqrt(mu/a^3) + ephem_struct.Delta_n;
    M_k = ephem_struct.M0 + n*prob.t_k;
    ecc = ephem_struct.e;
    tol = 1e-10; iter = 0; iter_max = 1000; % While loop control parameters:
    val_diff = 1; E_k = 1; % Initial Guess
    while abs(val_diff) > tol
        if iter > iter_max
            break
        end
        val_diff = (M_k-E_k+ecc*sin(E_k))/(-1+ecc*cos(E_k));
        E_k = E_k - val_diff;
        iter = iter+1;
    end
    sin_vk = sqrt(1-ecc^2)*sin(E_k)/(1-ecc*cos(E_k));
    cos_vk = (cos(E_k) - ecc)/(1-ecc*cos(E_k));
    v_k = atan2(sin_vk, cos_vk);
    phi_k_0 = v_k + ephem_struct.omega;
    phi_k = phi_k_0;
    d_phi_k = 1; tol = 1e-5; iter = 1; iter_max = 1000; % While loop control parameters
    while abs(d_phi_k) > tol
        if iter > iter_max
            break
        end
        d_phi_k = ephem_struct.Cus*sin(2*phi_k) + ephem_struct.Cuc*cos(2*phi_k);
        phi_k = phi_k_0 + d_phi_k;
        iter = iter+1;
    end
    u_k = phi_k;
    del_r_k = ephem_struct.Crs*sin(2*phi_k) + ephem_struct.Crc*cos(2*phi_k);
    del_i_k = ephem_struct.Cis*sin(2*phi_k) + ephem_struct.Cic*cos(2*phi_k);
    t = prob.t_k + ephem_struct.Toe;
    Omega_k = ephem_struct.OMEGA - O_dot_E*t + ephem_struct.OMEGA_DOT*prob.t_k;
    r_k = a*(1-ecc*cos(E_k)) + del_r_k;
    i_k = ephem_struct.i0 + ephem_struct.IDOT*prob.t_k + del_i_k;
    x_p = r_k*cos(u_k);
    y_p = r_k*sin(u_k);
    x_ecef = x_p*cos(Omega_k) - y_p*cos(i_k)*sin(Omega_k);
    y_ecef = x_p*sin(Omega_k) + y_p*cos(i_k)*cos(Omega_k);
    z_ecef = y_p*sin(i_k);
    sat_pos_calc = [x_ecef; y_ecef; z_ecef];
end

function prob_out = time_sort(prob_in)
    prob_out = struct();
    prob_out(1).x0 = []; % Receiver position
    prob_out(1).bu = []; % Receiver clock bias
    prob_out(1).G = []; % Satellite geometry matrix
    
    
    t_elapsed = [prob_in.TimeNanos];
    N_data = length(t_elapsed); % Total data points
    t_prev = t_elapsed(1); % Initial elapsed time
    ind = 1;
    ind_t = 1;
    ind_sat = 1;
    while ind < N_data
        while (t_elapsed(ind) - t_prev == 0)
            prob_out(ind_t).utcTimeMillis(ind_sat) = prob_in(ind).utcTimeMillis;
            prob_out(ind_t).TimeNanos(ind_sat) = prob_in(ind).TimeNanos;
            prob_out(ind_t).FullBiasNanos(ind_sat) = prob_in(ind).FullBiasNanos;
            prob_out(ind_t).BiasNanos(ind_sat) = prob_in(ind).BiasNanos;
            
            prob_out(ind_t).BiasUncertaintyNanos(ind_sat) = prob_in(ind).BiasUncertaintyNanos;
%             prob_out(ind_t).DriftNanosPerSecond(ind_sat) = prob_in(ind).DriftNanosPerSecond;
%             prob_out(ind_t).DriftUncertaintyNanosPerSecond(ind_sat) = prob_in(ind).DriftUncertaintyNanosPerSecond;
            prob_out(ind_t).HardwareClockDiscontinuityCount(ind_sat) = prob_in(ind).HardwareClockDiscontinuityCount;
            
            prob_out(ind_t).Svid(ind_sat) = prob_in(ind).Svid;
            prob_out(ind_t).TimeOffsetNanos(ind_sat) = prob_in(ind).TimeOffsetNanos;
            prob_out(ind_t).State(ind_sat) = prob_in(ind).State;
            prob_out(ind_t).ReceivedSvTimeNanos(ind_sat) = prob_in(ind).ReceivedSvTimeNanos;
            
            prob_out(ind_t).ReceivedSvTimeUncertaintyNanos(ind_sat) = prob_in(ind).ReceivedSvTimeUncertaintyNanos;
            prob_out(ind_t).Cn0DbHz(ind_sat) = prob_in(ind).Cn0DbHz;
            prob_out(ind_t).PseudorangeRateMetersPerSecond(ind_sat) = prob_in(ind).PseudorangeRateMetersPerSecond;
            prob_out(ind_t).PseudorangeRateUncertaintyMetersPerSecond(ind_sat) = prob_in(ind).PseudorangeRateUncertaintyMetersPerSecond;
            
            prob_out(ind_t).AccumulatedDeltaRangeState(ind_sat) = prob_in(ind).AccumulatedDeltaRangeState;
            prob_out(ind_t).AccumulatedDeltaRangeMeters(ind_sat) = prob_in(ind).AccumulatedDeltaRangeMeters;
            prob_out(ind_t).AccumulatedDeltaRangeUncertaintyMeters(ind_sat) = prob_in(ind).AccumulatedDeltaRangeUncertaintyMeters;
            prob_out(ind_t).CarrierFrequencyHz(ind_sat) = prob_in(ind).CarrierFrequencyHz;
            
            prob_out(ind_t).CarrierCycles(ind_sat) = prob_in(ind).CarrierCycles;
            prob_out(ind_t).MultipathIndicator(ind_sat) = prob_in(ind).MultipathIndicator;
%             prob_out(ind_t).AgcDb(ind_sat) = prob_in(ind).AgcDb;
            prob_out(ind_t).allRxMillis(ind_sat) = prob_in(ind).allRxMillis;
            
            prob_out(ind_t).sat_pos_calc(:,ind_sat) = prob_in(ind).sat_pos_calc;
            prob_out(ind_t).t_k(ind_sat) = prob_in(ind).t_k;
            prob_out(ind_t).B(ind_sat) = prob_in(ind).B;
            
            prob_out(ind_t).rho_meas(ind_sat) = rho_meas_calc(prob_in(ind));
            ind_sat = ind_sat + 1;
            ind = ind + 1;
            if ind > N_data, return; end
        end
        t_prev = t_elapsed(ind);
        ind_t = ind_t + 1;
        ind_sat = 1;
    end

end

function rho_m = rho_meas_calc(gnss_struct)
    week_s = 604800; c = 299792458;
    N_w = floor(-gnss_struct.FullBiasNanos/(week_s*1e9));
    t_Rx_hw = gnss_struct.TimeNanos + gnss_struct.TimeOffsetNanos;
    b_hw = gnss_struct.FullBiasNanos + gnss_struct.BiasNanos;
    t_Rx_GPS = t_Rx_hw-b_hw;
    t_Rx_w = t_Rx_GPS - N_w*week_s*1e9;
    rho_ns = t_Rx_w - gnss_struct.ReceivedSvTimeNanos;
    rho_m = rho_ns*c/1e9;
end

function G = G_mat(sat_pos, rx_pos, N)
    sat_rx = zeros(3,N);
    for ind = 1:N
        sat_pos_i = sat_pos(:,ind);
        sat_rx(:,ind) = -(sat_pos_i-rx_pos)/norm(sat_pos_i-rx_pos);  
    end
    G = [sat_rx' ones(N,1)];
end

function rho = compute_rho(sat_pos, rx_pos, b_u, B, N)
    rho = zeros(N,1);    
    for ind = 1:N
        sat_pos_i = sat_pos(:,ind);
        rho(ind,1) = norm(sat_pos_i-rx_pos)+b_u-B(ind);
    end
end

function prob = gyro_accel_mag_sort(accel, gyro, mag, prob)
    t_elapsed = [prob.utcTimeMillis];
    N_accel = length(accel);
    t_prev = t_elapsed(1); % Initial elapsed time
    ind = 1;
    ind_t = 1;
    ind_sat = 1;
    while ind < N_accel
        while (t_elapsed(ind) - t_prev == 0)
            prob(ind_t).a(ind_sat,:) = accel(ind,:);
            ind_sat = ind_sat + 1;
            ind = ind + 1;
            if ind > N_accel, break; end
        end
        t_prev = t_elapsed(ind);
        ind_t = ind_t + 1;
        ind_sat = 1;
    end
    
    N_gyro = length(gyro);
    t_prev = t_elapsed(1); % Initial elapsed time
    ind = 1;
    ind_t = 1;
    ind_sat = 1;
    while ind < N_gyro
        while (t_elapsed(ind) - t_prev == 0)
            prob(ind_t).w(ind_sat,:) = gyro(ind,:);
            ind_sat = ind_sat + 1;
            ind = ind + 1;
            if ind > N_gyro || ind > length(t_elapsed), break; end
        end
        if ind > N_gyro || ind > length(t_elapsed), break; end
        t_prev = t_elapsed(ind);
        ind_t = ind_t + 1;
        ind_sat = 1;
    end
    
    N_mag = length(mag);
    t_prev = t_elapsed(1); % Initial elapsed time
    ind = 1;
    ind_t = 1;
    ind_sat = 1;
    while ind < N_mag
        while (t_elapsed(ind) - t_prev == 0)
            prob(ind_t).mag(ind_sat,:) = mag(ind,:);
            ind_sat = ind_sat + 1;
            ind = ind + 1;
            if ind > N_mag || ind > length(t_elapsed), break; end
        end
        if ind > N_mag || ind > length(t_elapsed), break; end
        t_prev = t_elapsed(ind);
        ind_t = ind_t + 1;
        ind_sat = 1;
    end
end

function prob_sort = accel_calib(prob_sort)
    len_prob = length(prob_sort);
    for ind = 1:len_prob
        a_mat_uncal = prob_sort(ind).a;
        if isempty(a_mat_uncal)
            break;
        else
            a_mat_cal = [a_mat_uncal(:,3) - a_mat_uncal(:,6), ...
                a_mat_uncal(:,4) - a_mat_uncal(:,7),...
                a_mat_uncal(:,5) - a_mat_uncal(:,7)];
           prob_sort(ind).a_cal = a_mat_cal;
        end
    end
end

function prob_sort = gyro_calib(prob_sort)
    len_prob = length(prob_sort);
    for ind = 1:len_prob
        g_mat_uncal = prob_sort(ind).w;
        if isempty(g_mat_uncal)
            break;
        else
            a_mat_cal = [g_mat_uncal(:,3) - g_mat_uncal(:,6), ...
                g_mat_uncal(:,4) - g_mat_uncal(:,7),...
                g_mat_uncal(:,5) - g_mat_uncal(:,7)];
           prob_sort(ind).w_cal = a_mat_cal;
        end
    end
end

function prob_sort = mag_calib(prob_sort)
    len_prob = length(prob_sort);
    for ind = 1:len_prob
        m_mat_uncal = prob_sort(ind).mag;
        if isempty(m_mat_uncal)
            break;
        else
            a_mat_cal = [m_mat_uncal(:,3) - m_mat_uncal(:,6), ...
                m_mat_uncal(:,4) - m_mat_uncal(:,7),...
                m_mat_uncal(:,5) - m_mat_uncal(:,7)];
           prob_sort(ind).mag_cal = a_mat_cal;
        end
    end
end

function prob = NR_GPS(prob)
    tol = 1e-6; % Tolerance
    max_iter = 1000; % Maximum allowed iterations if non-convergent
    d_x = 1; iter = 0; % Initialize update and loop counter variables
    while norm(d_x)>tol
        iter = iter + 1;
        if iter > max_iter
            break
        end
        prob.G = G_mat(prob.sat_pos_unique, prob.x0, size(prob.sat_pos_unique,2));
        prob.rho_theo = compute_rho(prob.sat_pos_calc, prob.x0, prob.bu, prob.B, size(prob.sat_pos_calc,2));
        prob = rho_iono_free(prob, 'theo');
        d_rho = prob.rho_meas_IF'-prob.rho_theo_IF';
        d_x = (prob.G'*prob.G)\prob.G'*d_rho;
        prob.x0 = prob.x0 + d_x(1:3);
        prob.bu = prob.bu + d_x(end);
    end
end

function prob = rho_iono_free(prob,str)
    freq = prob.CarrierFrequencyHz;
    Svid = prob.Svid;
    if strcmp(str, 'meas')
        rho = prob.rho_meas;
    elseif strcmp(str, 'theo')
        rho = prob.rho_theo;
    else
        error('Incorrect argument')
    end
    unique_Svid = unique(Svid,'stable');
    unique_Svid_len = length(unique_Svid);
    rho_IF = zeros(1,unique_Svid_len);

    for Svid_ind = 1:unique_Svid_len
        ind_uni = find(Svid == unique_Svid(Svid_ind));
        if length(ind_uni) > 2
            error('More than two repeated measurements')
        elseif length(ind_uni) < 2
            rho_IF(Svid_ind) = rho(ind_uni);
        else
            F = freq(ind_uni);
            R = rho(ind_uni);
            rho_IF(Svid_ind) = F(1)^2*R(1)/(F(1)^2-F(2)^2) - F(2)^2*R(2)/(F(1)^2-F(2)^2);
        end      
    end
    if strcmp(str, 'meas')
        prob.rho_meas_IF = rho_IF;
        prob.Svid_unique = unique_Svid;
    elseif strcmp(str, 'theo')
        prob.rho_theo_IF = rho_IF;
    end
end

function prob = sat_pos_reduction(prob)
    unique_Svid_len = length(prob.Svid_unique);
    pos = zeros(3, unique_Svid_len);
    for Svid_ind = 1:unique_Svid_len
        ind_uni = find([prob.Svid] == prob.Svid_unique(Svid_ind));
        if length(ind_uni) > 2
            error('More than two repeated measurements')
        elseif length(ind_uni) < 2
            pos(:,Svid_ind) = prob.sat_pos_calc(:,ind_uni);
        else % Take the average of the two positions
            pos(:,Svid_ind) = mean(prob.sat_pos_calc(:,ind_uni), 2);
        end      
    end
    [prob.sat_pos_unique] = pos;
end
