%{
Filename: iEKF_solution.m
by, Shivam Soni - 06/05/2021
1) EKF - EKS Implementation
2) Run after "kalman_filtering.m"
%}
%% iEKF
clc;
% R and Q parameters:
r = .1;
q = 30;
% Initial State:
mu = [prob_sort(1).x0; 0; 0; 0; 0; 0; 0; 0];
clear prob_tight
tic
prob_tight = tightly_coupled(prob_sort, q, r, mu);
t = toc;
t/length(prob_tight)
%%
prob_tight(499).x0_kf_s = prob_tight(498).x0_kf_s;
prob_tight(500).x0_kf_s = prob_tight(498).x0_kf_s;
prob_tight(501).x0_kf_s = prob_tight(503).x0_kf_s;
prob_tight(502).x0_kf_s = prob_tight(503).x0_kf_s;

prob_tight(966).x0_kf_s = prob_tight(965).x0_kf_s;
prob_tight(967).x0_kf_s = prob_tight(968).x0_kf_s;

%% Plot LLA
x_lla = zeros(length(prob_tight),3); cmat = parula(length(prob_tight));
for ind=1:length(prob_tight)
    x_lla(ind,:) = ecef2lla(prob_tight(ind).x0_kf_s','WGS84'); 
end
figure;
for i=1:length(prob_tight)
    geoplot(x_lla(i,1), x_lla(i,2),'.','Color',cmat(i,:));hold on;
end
title('Receiver position in LLA frame using WGS84 datum');
colorbar('southoutside','TickLabelInterpreter','latex','FontSize',24,...
    'TicksMode','manual','Ticks',[0, 1], 'TickLabels',{'$t = 0$', '$t = t_{end}$'})
%% Save LLA
disp('Saving...')
writematrix(x_lla, 'LLA_EKS.csv')
disp('Saved!')
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
    
    T_tot = length(prob_sort); T_tot = 1295;
    for ind_t = 1:T_tot
        if isempty(prob_sort(ind_t).sat_pos_calc) || isempty(prob_sort(ind_t).a_cal) % No IMU or GNSS data
            % Use WLS previous position
            prob_out(ind_t).x0_kf = prob_sort(ind_t).x0;
            prob_out(ind_t).bu_kf = prob_sort(ind_t).bu;
            % Use previous timestep's values for the rest of the state variables
            prob_out(ind_t).v0_kf = prob_out(ind_t-1).v0_kf; 
            prob_out(ind_t).ang_kf = prob_out(ind_t-1).ang_kf;
            continue
        end
        
        if ind_t == T_tot
            dt = abs(prob_sort(ind_t).utcTimeMillis(1) - prob_sort(ind_t-1).utcTimeMillis(1))/1e6;
        else
            dt = abs(prob_sort(ind_t+1).utcTimeMillis(1) - prob_sort(ind_t).utcTimeMillis(1))/1e6;
        end

        % Control Matrices:
        A = diag(ones(10,1)) + diag([dt*ones(3,1); 0; 0; 0], 4);
        B = [zeros(4,6); diag(dt*ones(6,1))];

        % Measurement noise
        meas_len = length(prob_sort(ind_t).rho_meas);
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
        [mu, P] = ekf_forward(A, B, Q, R, P, mu, z, u_ecef', prob_sort(ind_t));

%         [mu, P] = kalman_filter_ekf(A, B, H, R, Q, P, mu, z, u_ecef', h);
        prob_out(ind_t).S = P;
        prob_out(ind_t).mu = mu;
        % Save updated states:
        if ind_t == 1
            prob_out(ind_t).x0_kf = mu(1:3);
            prob_out(ind_t).bu_kf = mu(4);
            prob_out(ind_t).v0_kf = mu(5:7);
            prob_out(ind_t).ang_kf = mu(8:10);
        elseif norm(prob_out(ind_t-1).x0_kf - mu(1:3)) > 1000 % 1km
            sprintf('AAAA')
            prob_out(ind_t).x0_kf = prob_sort(ind_t).x0; % Use WLS previous position
            prob_out(ind_t).bu_kf = prob_sort(ind_t).bu; 
            prob_out(ind_t).v0_kf = prob_out(ind_t-1).v0_kf; % Use previous timestep's values for the rest of the state variables
            prob_out(ind_t).ang_kf = prob_out(ind_t-1).ang_kf;
        else
            prob_out(ind_t).x0_kf = mu(1:3);
            prob_out(ind_t).bu_kf = mu(4);
            prob_out(ind_t).v0_kf = mu(5:7);
            prob_out(ind_t).ang_kf = mu(8:10);
        end
    end
    % Backward Pass:
    % Initialize:
    mu_s = prob_out(T_tot).mu;
    S_s = prob_out(T_tot).S;
    for ind_t = T_tot:-1:1
        if isempty(prob_sort(ind_t).sat_pos_calc) || isempty(prob_sort(ind_t).a_cal) % No IMU or GNSS data
            % Use WLS previous position
            prob_out(ind_t).x0_kf = prob_sort(ind_t).x0;
            prob_out(ind_t).bu_kf = prob_sort(ind_t).bu;
            % Use previous timestep's values for the rest of the state variables
            prob_out(ind_t).v0_kf = prob_out(ind_t-1).v0_kf; 
            prob_out(ind_t).ang_kf = prob_out(ind_t-1).ang_kf;
            continue
        end
        if ind_t == T_tot
            dt = abs(prob_sort(ind_t).utcTimeMillis(1) - prob_sort(ind_t-1).utcTimeMillis(1))/1e6;
        else
            dt = abs(prob_sort(ind_t+1).utcTimeMillis(1) - prob_sort(ind_t).utcTimeMillis(1))/1e6;
        end
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
        % Control Matrices:
        A = diag(ones(10,1)) + diag([dt*ones(3,1); 0; 0; 0], 4);
        B = [zeros(4,6); diag(dt*ones(6,1))];
        
        % Backward Pass:
        [mu_s, S_s] = eks_backward(A, B, prob_out(ind_t).mu, prob_out(ind_t).S,...
            mu_s, S_s, u_ecef', Q);
        prob_out(ind_t).mu_s = mu_s;
        prob_out(ind_t).S_s = S_s;
        if ind_t == 1
            prob_out(ind_t).x0_kf_s = mu_s(1:3);
            prob_out(ind_t).bu_kf_s = mu_s(4);
            prob_out(ind_t).v0_kf_s = mu_s(5:7);
            prob_out(ind_t).ang_kf_s = mu_s(8:10);
        elseif norm(prob_out(ind_t-1).x0_kf - mu_s(1:3)) > 1000 % 1km
            sprintf('AAAA')
            prob_out(ind_t).x0_kf_s = prob_sort(ind_t).x0; % Use WLS previous position
            prob_out(ind_t).bu_kf_s = prob_sort(ind_t).bu; 
            prob_out(ind_t).v0_kf_s = prob_out(ind_t-1).mu(5:7); % Use previous timestep's values for the rest of the state variables
            prob_out(ind_t).ang_kf_s = prob_out(ind_t-1).mu(8:10);
        else
            prob_out(ind_t).x0_kf_s = mu_s(1:3);
            prob_out(ind_t).bu_kf_s = mu_s(4);
            prob_out(ind_t).v0_kf_s = mu_s(5:7);
            prob_out(ind_t).ang_kf_s = mu_s(8:10);
        end
    end
end


function [mu_t_t, S_t_t] = ekf_forward(A, B, Q, R, S_tm_tm, mu_tm_tm, y, u, prob_sort)
    % Predict:
    mu_t_tm = A*mu_tm_tm + B*u; % Mean
    S_t_tm = A*S_tm_tm*A' + Q; % Covariance
    
    % Measurement Jacobian:
    C = meas_J(mu_t_tm, prob_sort);
    
    % Update:
    K_t = S_t_tm*C'/(C*S_t_tm*C'+R);
    mu_t_t = mu_t_tm + K_t*(y - meas_fn(mu_t_tm, prob_sort));
    S_t_t = S_t_tm - K_t*C*S_t_tm;
end

function [mu_t_T, S_t_T] = eks_backward(A, B, mu_t_t, S_t_t, mu_tp_T, S_tp_T, u, Q)
    % Predict:
    mu_tp_t = A*mu_t_t + B*u; % Mean
    S_tp_t = A*S_t_t*A' + Q; % Covariance
    
    % Smooth:
    K_t_s = S_t_t*A'/S_tp_t;
    mu_t_T = mu_t_t + K_t_s*(mu_tp_T - mu_tp_t);
    S_t_T = S_t_t + K_t_s*(S_tp_T - S_tp_t)*K_t_s';
end


function C = meas_J(mu, prob_sort)
    meas_len = length(prob_sort.rho_meas);
    C = zeros(meas_len, 10);
    for ind_sat = 1:meas_len
        X = prob_sort.sat_pos_calc(1,ind_sat) - mu(1);
        Y = prob_sort.sat_pos_calc(2,ind_sat) - mu(2);
        Z = prob_sort.sat_pos_calc(3,ind_sat) - mu(3);
        eta = sqrt(X^2 + Y^2 + Z^2);
        C(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
    end
end

function h = meas_fn(mu, prob_sort)
    meas_len = length(prob_sort.rho_meas);
    h = zeros(meas_len, 1);
    for ind_sat = 1:meas_len
        X = prob_sort.sat_pos_calc(1,ind_sat) - mu(1);
        Y = prob_sort.sat_pos_calc(2,ind_sat) - mu(2);
        Z = prob_sort.sat_pos_calc(3,ind_sat) - mu(3);
        eta = sqrt(X^2 + Y^2 + Z^2);
        h(ind_sat,1) = eta + mu(4) - prob_sort.B(ind_sat);
    end
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
