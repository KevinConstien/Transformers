clc
close all
clear all
%IF 1000, Average 25
filename(1) = "2019-03-25_20_32_50AntennaDistance21cm.csv";
filename(2) = "2019-03-25_20_13_58AntennaDistance24cm.csv";
filename(3) = "2019-03-25_20_04_58AntennaDistance27cm.csv";
filename(4) = "2019-03-25_19_56_06AntennaDistance30cm.csv";
filename(5) = "2019-03-25_18_22_30AntennaDistance33cm.csv";
colors = ["r-","k-","g-","b-","c-"];
labels = ["33cm","30cm","27cm","24cm","21cm"];
% otherLabels = ["85\lambda","77\lambda","69\lambda","62\lambda","54\lambda"];

for i = 1:length(filename)
    coTable = ImportRwave(filename(i));

    freq = 7.70e10;
    rows = (coTable.FreqHz == freq);
    r = linspace(-45,45,91);


    S21 = coTable.S21DB(rows);
    S12 = coTable.S12DB(rows);
    S21magPhaseComparison=coTable.S21DB(rows);
    S21anglePhaseComparison = coTable.S21DEG(rows);
    % delta = (abs(max(S21))-abs(max(S12)));
    % S21 = S21 + delta;
    maxs21 = max(S21);




    
    plot(r,S21,colors(i),'linewidth',3); hold on
    
    %plot([-60:1:60],S21angle,'r-','linewidth',3);
    %plot([-60:1:60],newabsorbers1degree,'k', 'linewidth',3);
    %plot([-60:0.5:60],oldabsorbers,'b-','linewidth',3);
    %plot(coTable.y(rows),S12,'r-','linewidth',3);

end

grid on 
grid minor
legend(labels(:));
set(gca,'FontSize', 18);
title('FF Horizontal Cut');
ylabel('Magnitude (dB)','fontsize',18);
xlabel('Angle (degrees)','fontsize',18);
freqq = num2str(freq/1e9);
maximum = num2str(maxs21,3);
    
%     text(-45,-30,{['Freq: ', freqq,' GHz'],['S21(max): ', maximum,' dB']},'fontsize',12 );
