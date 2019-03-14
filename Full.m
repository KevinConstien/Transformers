clc
close all
clear all


%Importing H data
% filename = '2019-03-08_00_40_29test2NF.csv';
filename = '2019-03-13_20_30_25testNFAver25IF3kFoamOn.csv';
coTable = ImportRwave(filename);

save('coTable.mat','coTable');

% load 'coTable.mat';

%Importing V data
% filename = '2019-03-08_00_40_29test2NF.csv';
cxTable = coTable;
% cxTable = ImportRwave(filename);

%Whatever freauency you want. I supposed you could just loop for all freqs
% freq = 5.3450e9;
freq = 75000000000;
rows = (coTable.FreqHz == freq);

%Reformatting data for transformation
co(:,1) = coTable.x(rows);
co(:,2) = coTable.y(rows);
co(:,3) = coTable.S21DB(rows);
co(:,4) = coTable.S21DEG(rows);

cross(:,1) = cxTable.x(rows);
cross(:,2) = cxTable.y(rows);
cross(:,3) = cxTable.S21DB(rows);
cross(:,4) = cxTable.S21DEG(rows);

%Transformation Parameters
FFparams.pei = true;
FFparams.freq = freq;
FFparams.nbr_samples = 361;
FFparams.H_angle_range = 120;
FFparams.V_angle_range = 120;
[FF] = nf2ffFunction(FFparams,co,cross);
% 
% %Plotting
FF.coordinate = 'L3';
% kplot(FF);











