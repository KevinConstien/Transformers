function [nbr_cases] = mainLoop(indexPattern)
%% load data
% filename = '2019-04-11_23_14_57NFScanAt15mmWDiffAbsor.csv';
% coTable = ImportRwave(filename);
% save('coTable.mat','coTable');

load('coTable.mat')
cxTable = coTable;

nbr_freqs = 10;
truncation = 2;
%%
%Whatever freauency you want. I supposed you could just loop for all freqs
% freq = 5.3450e9;
[nbr_samples_total,~] = size(coTable);  

% freq = 77000000000;
freq = 77e9;
rows = (coTable.FreqHz == freq);

[nbr_rows_per_freq, ~]=size(coTable.x(rows));
rows = find(rows);

% indexPattern =11;


%calculate how many patterns were taken at the same position
nbr_cases = sum(coTable.x(rows) == coTable.x(1) & coTable.y(rows) == coTable.y(1));
nbr_PointsInX     = sum(coTable.x(rows) == coTable.x(1))/nbr_cases;
nbr_PointsInY     = sum(coTable.y(rows) == coTable.y(1))/nbr_cases;

%add a Case_ID column to the table
ID = [1:nbr_cases];
Case_ID =repmat(ID,nbr_freqs,1); %per frequencies
Case_ID_col = repmat(Case_ID(:),nbr_PointsInY*nbr_PointsInX,1);
size(Case_ID_col);
coTable = [coTable table(Case_ID_col)];
coTable = sortrows(coTable,[10 11],'descend'); %this is because the scan raster

%trim data for a case
clearvars 'rows'
rows = find(coTable.FreqHz == freq & coTable.Case_ID_col==indexPattern);
[nbr_rows, ~]=size(coTable.x(rows));

%grid
X = -reshape(coTable.x(rows),sqrt(nbr_rows),sqrt(nbr_rows));%repmat(grid,sqrt(nbr_rows),1);
Y = -reshape(coTable.y(rows),sqrt(nbr_rows),sqrt(nbr_rows));%repmat(grid',1,sqrt(nbr_rows));%

%apply window
sampleCoordinates = X(1,:);
index2sample = find(...
    X >= sampleCoordinates(truncation+1) &  X <= sampleCoordinates(nbr_PointsInX-truncation) &...
    Y >= sampleCoordinates(truncation+1) &  Y <= sampleCoordinates(nbr_PointsInX-truncation));

% %plot
S21db = reshape(coTable.S21DB(rows),sqrt(nbr_rows),sqrt(nbr_rows));
S21deg = reshape(coTable.S21DEG(rows),sqrt(nbr_rows),sqrt(nbr_rows));
% surf(X,Y,S21db)
% view(2)
% shading interp
% colormap jet
% colorbar
% saveas(gcf,'Amplitude','png')
% 
% figure
% surf(X,Y,S21deg)
% view(2)
% shading interp
% colormap jet
% colorbar
% saveas(gcf,'Phase','png')

%Reformatting data for transformation

co(:,1) = X(index2sample);
co(:,2) = Y(index2sample);
co(:,3) = S21db(index2sample);
co(:,4) = S21deg(index2sample);

cross = co;
cross(:,3) = cross(:,3)./cross(:,3)*-300;
% cross(:,1) = cxTable.x(rows);
% cross(:,2) = cxTable.y(rows);
% cross(:,3) = cxTable.S12DB(rows);%./abs(cxTable.S12DB(rows))*-200;
% cross(:,4) = cxTable.S12DEG(rows);

%%
%Transformation Parameters
FFparams.pei = true;
FFparams.freq = freq;
FFparams.nbr_samples = 361;
FFparams.H_angle_range = 100;
FFparams.V_angle_range = 100;
[FF] = nf2ffFunction(FFparams,co,cross);

%Plotting
FF.coordinate = 'L3';
kplotGIF(FF);
end











