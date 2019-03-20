
%% Import data from text file.
% Script for importing data from the following text file:
%
%    G:\My Drive\ARRC\NF2FF\nf2ff\KevinCode\InProgress\sampleH.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.
 
% Auto-generated by MATLAB on 2019/02/11 10:51:29
% clc
% close all
% clear all
 
 
 
 
%% Initialize variables.
% filename = '2019-01-22_21_59_46_00_H_UniformUniform.csv';
function [data] = tempThing(filename)
fileID = fopen(filename,'r');
all = ((textscan(fileID,'%s%s%s%s%s%s%s%s%s%[^\n\r]','Delimiter',',')));
fclose(fileID);
firstcol = string(all{:,1});
start = firstcol(1);
AOlength = 64;
 
%Detecting the length for each chunck of frequencies
for i = 2:length(firstcol)
    if firstcol(i) == start
        break
    end
end
 
fileID = fopen(filename,'r');
headerlines = 7;
delimiter = ',';
chunk = i-headerlines-1;
numCols = 12+64;
 

while ~feof(fileID) 
%     format = ',%d %s %s  %f  %f %f %f %f %s %s';
    tic
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    header = textscan(fileID,formatSpec,headerlines,'Delimiter',',','Headerlines',0);
    flathead = [header{:}];
    if feof(fileID)
        break
    end
     
    ypos = string(flathead(4,1));
    xpos = string(flathead(4,2));
    rpos = string(flathead(4,3));
    x(1:chunk,1) = xpos;y(1:chunk,1) = ypos; r(1:chunk,1) = rpos;
     
    A0vals = zeros(chunk,AOlength);
    for i = 1:AOlength
        A0vals(1:chunk,i) = string(flathead(6,i));
    end
 
 
    formatSpec = '%s%s%s%s%s%s%s%s%s%[^\n\r]';
    tempArray = textscan(fileID, formatSpec,chunk, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,0, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    tempArray{1,10} = string(x); tempArray{1,11} = string(y); tempArray{1,12} = string(r);
    %RODRIGO - add AO chunk to tempArray
    lastIndex = 12;%index for r
    for i=1:AOlength
        tempArray{1,i+lastIndex} = string(A0vals(:,i));
    end
    
    if exist('dataArray','var') == 0
        dataArray = tempArray;
    else
        dataArray = cat(1,dataArray,tempArray);
    end
    clear tempArray

    if toc > .02
        if exist('storage','var') == 0
            storage = dataArray;
        else
            storage = cat(1,storage,dataArray);
        end
        clear dataArray;
    end
end

if exist('storage','var') == 0;
    myflatcellarray = [dataArray{:}];
else
    myflatcellarray=[storage{:}];
end

sample = (reshape(myflatcellarray,[],numCols));
 
 
%% Close the text file.
fclose(fileID);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KEVIN'S MUCH QUICKER SOLUTION. IF THE DATA IS NOT BEING 
%REPRESENTED CORRECTLY, COMMENT OUT AND USE THE SECTION
%BELOW INSTEAD.
doubleCols = 1:numCols;
for col = 1:numCols
    for row = 1:length(sample)
        if sum(col == doubleCols)
            raw{row,col} = double(sample(row,col));
        else
            raw{row,col} = sample(row,col);
        end
    end    
end
 
 
 
%% Create output variable
data = table;
data.FreqHz = cell2mat(raw(:, 1));
data.S11DB = cell2mat(raw(:, 2));
data.S11DEG = cell2mat(raw(:, 3));
data.S21DB = cell2mat(raw(:, 4));
data.S21DEG = cell2mat(raw(:, 5));
data.S12DB = cell2mat(raw(:, 6));
data.S12DEG = cell2mat(raw(:, 7));
data.S22DB = cell2mat(raw(:, 8));
data.S22DEG = cell2mat(raw(:, 9));
data.x = cell2mat(raw(:, 10));
data.y = cell2mat(raw(:, 11));
data.r = cell2mat(raw(:, 12));
 
%Rodrigo - Analog Outputs
TAO = cell2table(raw(:,(lastIndex+1):end));
for i=1:AOlength %rename cols
    if i <= 32
        colName = strcat('B2','AO',num2str(i-1));
        TAO.Properties.VariableNames(i) = {colName};
    elseif i <= 64
        colName = strcat('B1','AO',num2str(i-32));
        TAO.Properties.VariableNames(i) = {colName};
    end
end
data = [data TAO];
 
%Rename
save('data.mat','data');
end
 
%% Clear temporary variables
% clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;