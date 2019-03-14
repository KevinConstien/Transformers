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
    S21 = coTable.S21DB(rows);
    S21angle = coTable.S21DEG(rows);

    figure()
    plot(S21);

    title("FF Horizontal Cut, Freq = " + freq/1e9+"GHz");
    ylabel('Magnitude (dB)');
    xlabel('Angle (\phi)');
end



