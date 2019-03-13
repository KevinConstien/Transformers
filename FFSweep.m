clc
close all
clear all

% filename = '2019-03-13_20_30_25testNFAver25IF3kFoamOn.csv';
coTable = ImportRwave(filename);

freq = 75000000000;
rows = (coTable.FreqHz == freq);


r = coTable.r(rows);
S21 = coTable.S21DB(rows);
S21angle = coTable.S21DEG(rows);

figure()
plot(r,S21);

title('FF Horizontal Cut');
xlabel('Magnitude (dB)');
ylabel('Angle (\phi)');


