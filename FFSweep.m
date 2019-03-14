clc
close all
clear all

filename = '2019-03-13_22_45_13testFFArrayIF3kAver25.csv';
coTable = ImportRwave(filename);

freq = 75000000000;


for i = 1:20:length(unique(coTable.FreqHz))
    freq = coTable.FreqHz(i);
    rows = (coTable.FreqHz == freq);
    r = coTable.y(rows);

    S12 = coTable.S12DB(rows);
    S21 = coTable.S21DB(rows);   

    figure();
    plot(r - r(end)/2,S12,'-o');hold on;
    plot(r - r(end)/2,S21,'-o');

    title("FF Horizontal Cut, Freq = " + freq/1e9+"GHz");
    legend('S_1_2','S_2_1');
    ylabel('Magnitude (dB)');
    xlabel('Angle (\phi)');
end



