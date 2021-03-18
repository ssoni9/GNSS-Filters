accel = table2array(readtable('accel.csv','ReadVariableNames',0));
gyro = table2array(readtable('gyro.csv','ReadVariableNames',0));
prob_sort(1).a = [];
prob_sort(1).w = [];
prob_sort = gyro_accel_sort(accel, gyro, prob_sort);

%% Kalman Filter:
% R and Q parameters:
r = 30;
q = 5;
% Initial State:
mu = [prob_sort(1).x0; 0; 0; 0; 0; 0; 0; 0];
clear prob_tight
prob_tight = tightly_coupled(prob_sort, q, r, mu);


%% Plot LLA
x_lla = zeros(length(prob_tight),3); cmat = parula(length(prob_tight));
for ind=1:length(prob_tight)
    x_lla(ind,:) = ecef2lla(prob_tight(ind).x0_kf','WGS84'); 
end
figure;
for i=1:length(prob_tight)
    geoplot(x_lla(i,1), x_lla(i,2),'.','Color',cmat(i,:));hold on;
end
title('Receiver position in LLA frame using WGS84 datum');
colorbar('southoutside','TickLabelInterpreter','latex','FontSize',24,...
    'TicksMode','manual','Ticks',[0, 1], 'TickLabels',{'$t = 0$', '$t = t_{end}$'})

%% Functions:

function prob_out = tightly_coupled(prob_sort, q, r, mu)
    % Output Problem:
    prob_out = struct();
    
    % Initialize Covariance: 
    P = diag(rand(10,1));
    
    % Process Noise:
    Q = q*diag(ones(10,1));
    
    % State:
    prob_out(1).x0_kf = mu(1:3);
    prob_out(1).bu_kf = mu(4);
    prob_out(1).v0_kf = mu(5:6);
    prob_out(1).ang_kf = mu(7:9);
    
    T_tot = length(prob_sort);
    for ind_t = 1:T_tot
        if isempty(prob_sort(ind_t).sat_pos_calc) || isempty(prob_sort(ind_t).a_cal)
            break
        end
        
        if ind_t == T_tot
            dt = abs(prob_sort(ind_t).utcTimeMillis(1) - prob_sort(ind_t-1).utcTimeMillis(1))/1e6;
        else
            dt = abs(prob_sort(ind_t+1).utcTimeMillis(1) - prob_sort(ind_t).utcTimeMillis(1))/1e6;
        end

        % Control Matrices:
        A = diag(ones(10,1)) + diag([dt*ones(3,1); 0; 0; 0], 4);
        B = [zeros(4,6); diag(dt*ones(6,1))];
        
        % Jacobian:
        meas_len = length(prob_sort(ind_t).rho_meas);
        H = zeros(meas_len, 10); h = zeros(meas_len, 1);
        for ind_sat = 1:meas_len
            if ind_t == 1
                X = prob_sort(ind_t).sat_pos_calc(1,ind_sat) - prob_out(ind_t).x0_kf(1);
                Y = prob_sort(ind_t).sat_pos_calc(2,ind_sat) - prob_out(ind_t).x0_kf(2);
                Z = prob_sort(ind_t).sat_pos_calc(3,ind_sat) - prob_out(ind_t).x0_kf(3);
            else
                X = prob_sort(ind_t).sat_pos_calc(1,ind_sat) - prob_out(ind_t-1).x0_kf(1);
                Y = prob_sort(ind_t).sat_pos_calc(2,ind_sat) - prob_out(ind_t-1).x0_kf(2);
                Z = prob_sort(ind_t).sat_pos_calc(3,ind_sat) - prob_out(ind_t-1).x0_kf(3);
            end    
            eta = sqrt(X^2 + Y^2 + Z^2);
            H(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
            if ind_t == 1
                h(ind_sat,1) = eta + prob_out(ind_t).bu_kf - prob_sort(ind_t).B(ind_sat);
            else
                h(ind_sat,1) = eta + prob_out(ind_t-1).bu_kf - prob_sort(ind_t).B(ind_sat);
            end
        end

        % Measurement noise
        R = r*diag(ones(meas_len,1));

        % Measurement:
        z = [prob_sort(ind_t).rho_meas]';

        % Disturbances or Acceleration Inputs:
        acc = prob_sort(ind_t).a_cal;
        gyr = prob_sort(ind_t).w_cal;
        mag = prob_sort(ind_t).mag_cal;
        [u_bod, mag] = kf_meas_vec(acc, gyr, mag);
        if ind_t == 1
            u_ecef = body2ecef(u_bod, mag, prob_out(ind_t).x0_kf);
        else
            u_ecef = body2ecef(u_bod, mag, prob_out(ind_t-1).x0_kf);
        end
        
        % Call Kalman Filter
        [mu, P] = kalman_filter_ekf(A, B, H, R, Q, P, mu, z, u_ecef', h);

        % Save updated states:
        prob_out(ind_t).x0_kf = mu(1:3);
        prob_out(ind_t).bu_kf = mu(4);
        prob_out(ind_t).v0_kf = mu(5:7);
        prob_out(ind_t).ang_kf = mu(8:10);
    end
end

function u_ecef = body2ecef(acc_gyr, mag, pos)
    acc = acc_gyr(1:3);
    gyr = acc_gyr(4:6);
    pos = reshape(pos, 1, 3);
    gvec_ned = [0, 0, 9.81];
    orientation = ecompass(acc, mag); % a_ned_d = acos(9.81/a_y)
    avec_ned = rotatepoint(orientation, acc);
    avec_ned_no_g = avec_ned - gvec_ned;
    lla = ecef2lla(pos, 'WGS84');
    R_ecef2ned = RotEcef2Ned(lla(1), lla(2));
    avec_ecef = R_ecef2ned'*avec_ned_no_g';
    wvec_ecef = R_ecef2ned'*gyr';
    u_ecef = [reshape(avec_ecef, 1, 3), reshape(wvec_ecef, 1, 3)];
end

function [mu_t_t, P_t_t] = kalman_filter_ekf(A, B, H, R, Q, P_tm_tm, mu_tm_tm, z_t, u_t, h_mu_t_tm)
    % Predict:
    mu_t_tm = A*mu_tm_tm + B*u_t;
    P_t_tm = A*P_tm_tm*A' + Q;
    
    % Update:
    y_t = z_t - h_mu_t_tm;
    K_t = P_t_tm*H'/(R + H*P_t_tm*H');
    mu_t_t = mu_t_tm + K_t*y_t;
    temp = K_t*H; temp = eye(size(temp)) - temp;
    P_t_t = temp*P_t_tm*temp' + K_t*R*K_t';
end

function [u_vec, m_sm] = kf_meas_vec(a, w, m)
    a_sm = smooth_mean(a);
    w_sm = smooth_mean(w);
    m_sm = smooth_mean(m);
    u_vec = [reshape(a_sm, 1, 3), reshape(w_sm, 1, 3)];
end


function vec_out = smooth_mean(mat)
    vec_out = mean(smoothdata(mat, 1, "gaussian", [4,4]), 1); 
end