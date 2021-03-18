%{
Filename: data_manip.m
by, Shivam Soni - 03/12/2021
1) Reads the raw csv file and outputs a "gnss_log.csv" table
2) Reads the hourly gps ephemeris file and outputs "ephem_file.csv" table
%}
%%
clearvars; clc; close all;

% Read the raw data:
Data_struct = ReadGnssLogger(pwd, 'Mi8_GnssLog.txt', SetDataFilter);
Data_table = struct2table(Data_struct);

% Save the raw data:
writetable(Data_table, 'gnss_log.csv');

% Read hourly ephemeris Rinex file:
Eph_struct = ReadRinexNav('hour1990.20n');

% Save the hourly ephemeris file:
Eph_table = struct2table(Eph_struct);
writetable(Eph_table, 'ephem_file.csv')
%% Read accelerations:
clearvars; clc; close all
accel_filename = 'Mi8_Accel.txt';
Data_cell = readcell(accel_filename);

len_data = length(Data_cell);
Data_out = struct(); a_ind = 1; g_ind = 1; m_ind = 1;
for ind = 1:len_data
    if strcmp(Data_cell(ind,1), 'Accel')
        Data_out.accel(a_ind,:) = cell2mat(Data_cell(ind,2:9));
        a_ind = a_ind + 1;
    elseif strcmp(Data_cell(ind,1), 'Gyro')
        Data_out.gyro(g_ind,:) = cell2mat(Data_cell(ind,2:9));
        g_ind = g_ind + 1;
    elseif strcmp(Data_cell(ind,1), 'Mag')
        Data_out.mag(m_ind,:) = cell2mat(Data_cell(ind,2:9));
        m_ind = m_ind + 1;
    else
        continue
    end
end
%% Save as csv:
writematrix(Data_out.accel, 'accel.csv');
writematrix(Data_out.gyro, 'gyro.csv');
writematrix(Data_out.mag, 'mag.csv');
