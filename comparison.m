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

[~,~,h_err] = Nmea2ErrorPlot('test_EKF.nmea','SPAN_Mi8_10Hz.nmea', pwd);

%% Plot states
angles = [prob_tight.ang_kf];

figure; hold on;
subplot(3, 1, 1); plot(angles(1,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\psi$ [deg]'); legend('$\psi$');

subplot(3, 1, 2); plot(angles(2,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\theta$ [deg]');  legend('$\theta$');

subplot(3, 1, 3); plot(angles(3,:), 'LineWidth', 2); grid minor; 
xlabel('\# measurement'), ylabel('$\phi$ [deg]');  legend('$\phi$');



