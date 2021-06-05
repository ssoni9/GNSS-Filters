%% 
clearvars; clc; close all; 
%% Read Ground Truth CSV Files:
gga = table2array(readtable('GT_gga.csv','ReadVariableNames',1));
rmc = table2array(readtable('GT_rmc.csv','ReadVariableNames',1));

%% 
%{
Conduct UTC [ms] to UTC [yyyy-mm-dd-hh-mm-ss] using the following online
resource: https://currentmillis.com/

Convert the LLA file generated using WLS or EKF using the following online
resource: https://www.nmeagen.org/
%}
%% Generate Plots using Google's functions:

[~,~,h_err] = Nmea2ErrorPlot('test_EKS.nmea','SPAN_Mi8_10Hz.nmea', pwd);
%% Revised Plot:
n = 45; % Truncate outliers:
temp = h_err.distM(h_err.distM < n);
dM = sort(temp);
figure; hold on; HH = cdfplot(dM); XD = HH.XData; YD = HH.YData;
clf; plot(XD, YD,'LineWidth', 2);
ts = sprintf('cdf (%d points)',length(dM));
title(ts)
xlabel('Horizontal error [m]')
ylabel('Fraction of distribution')
hold on

dM50 = median(dM);
plot(dM50,.5,'.r','MarkerSize', 16)
ts = ['$ 50\% $ ', num2str(dM50)];
text(dM50,.5,ts, 'VerticalAlignment','Top', 'Interpreter', 'latex', 'FontSize', 16)

dM67 = dM(round(length(dM)*.67));
plot(dM67,.67,'.r','MarkerSize', 16)
ts = ['$ 67\% $ ', num2str(dM50)];
text(dM67,.67,ts,'VerticalAlignment','Top', 'Interpreter', 'latex', 'FontSize', 16)

dM95 = dM(round(length(dM)*.95));
plot(dM95,.95,'.r','MarkerSize', 16)
ts = ['$ 95\% $ ', num2str(dM95)];
text(dM95,.95,ts,'VerticalAlignment','Top', 'Interpreter', 'latex', 'FontSize', 16)

dM100 = max(dM);
plot(dM100,1,'.r','MarkerSize', 16)
ts = ['$ 100\% $ ', num2str(dM100)];
text(dM100-3,1-0.01,ts,'VerticalAlignment','Top', 'Interpreter', 'latex', 'FontSize', 16)

grid minor
title('EKF-EKS: Horizontal Error CDF Plot')

%% Plot states
angles = [prob_tight.ang_kf];

figure; hold on;
subplot(3, 1, 1); plot(angles(1,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\psi$ [deg]'); legend('$\psi$');

subplot(3, 1, 2); plot(angles(2,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\theta$ [deg]');  legend('$\theta$');

subplot(3, 1, 3); plot(angles(3,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\phi$ [deg]');  legend('$\phi$');



