function [FFexport] = nf2ffFunction(FFparams,co,cross)

%If x,y coordinates not in terms of meters, then convert!
co(:,1:2) = co(:,1:2) /1e3;
cross(:,1:2) = cross(:,1:2)/1e3;

% Need to know the dimensions of the position data in order to
% Properly reshape. May need to change this method in the future
% Depending on how the imported data looks.
numbersx = length(unique(rmmissing(co(:,1))));
numbersy = length(unique(rmmissing(co(:,2))));
%Creating complex vectors
co_cmplx = db2mag(co(:,3)).*exp(1i*co(:,4)*pi/180);
cx_cmplx = db2mag(cross(:,3)).*exp(1i*cross(:,4)*pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Turning NF data into a 2D grid
%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(FFparams,'pei') && FFparams.pei == true
    gridx = reshape(co(:,1),[],numbersx);
    gridy = reshape(co(:,2),[],numbersx);
    
%     totalx = 15*28e-3;
%     x = linspace(-totalx/2,totalx/2,16);
%     [gridx,gridy] = meshgrid(x,x);

    gridco = reshape(co(:,3),[],numbersx);
    gridpco = reshape(co(:,4),[],numbersx);
    gridco_cmplx = reshape(co_cmplx,[],numbersx);

    gridcx = reshape(cross(:,3),[],numbersx);
    gridpcx = reshape(cross(:,4),[],numbersx);
    gridcx_cmplx = reshape(cx_cmplx,[],numbersx);
else 
    gridx = reshape(co(:,1),[],numbersy).';
    gridy = reshape(co(:,2),[],numbersy).';
    gridco = reshape(co(:,3),[],numbersy);gridco = gridco.';
    gridpco = reshape(co(:,4),[],numbersy);gridpco = gridpco.';
    gridco_cmplx = reshape(co_cmplx,[],numbersy);gridco_cmplx = gridco_cmplx.';
    gridcx = reshape(cross(:,3),[],numbersy);gridcx = gridcx.';
    gridpcx = reshape(cross(:,4),[],numbersy);gridpcx = gridpcx.';
    gridcx_cmplx = reshape(cx_cmplx,[],numbersy);gridcx_cmplx = gridcx_cmplx.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth = 2;

%%% Simple derivations
lambda = physconst('LightSpeed')/FFparams.freq;
[Nx, My] = size(gridx);
NI = 128; %nbr freqs in y
MI = 128; %nbr freqs in x
dx = abs(gridx(1,2)-gridx(1,1));
dy = abs(gridy(2,1)-gridy(1,1));


%%% Process NF2FF
%create grid of FF sample angles --------------
H_angles = linspace(-FFparams.H_angle_range/2,FFparams.H_angle_range/2,FFparams.nbr_samples);
V_angles = linspace(-FFparams.V_angle_range/2,FFparams.V_angle_range/2,FFparams.nbr_samples);

[AZ,EL] = meshgrid(H_angles,V_angles);
[U,V,W] = AzEl2UVW(AZ,EL);
[TH,PH] = UVW2ThPh(U,V,W);

%Compute FF --------------
K = 2*pi/lambda;
m = -MI/2:1:(MI/2-1);
n = -NI/2:1:(NI/2-1); %from balanis
[M,N] = meshgrid(m,n);
Kx = 2*K*M/MI;
Ky = 2*K*N/NI;
Kz = sqrt(K^2 - Kx.^2 - Ky.^2);


%DFT -----------
%deltas for trapezoidal approximation NARASIMHAM 1984

Fco = zeros(NI,MI);
Fcross = zeros(NI,MI);

for i=1:NI
    for j=1:MI
        Fco(i,j) = sum(sum(gridco_cmplx.*exp(1i*(Kx(i,j)*gridx+Ky(i,j)*gridy))*dx*dy));
        Fcross(i,j)  = sum(sum(gridcx_cmplx .*exp(1i*(Kx(i,j)*gridx+Ky(i,j)*gridy))*dx*dy));
    end
end

Fco_Interp = interp2(Kx, Ky, Fco,K*U,K*V,'spline');
Fcross_Interp = interp2(Kx, Ky, Fcross,K*U,K*V,'spline');

pol = 1;
if pol==1 %H,AZ,X pol
    F_AZ = Fco_Interp;
    F_EL = Fcross_Interp;
elseif pol==0 %V,EL or Y
    F_AZ = Fcross_Interp;
    F_EL = Fco_Interp;
end

%PROBE CORRECTION ---------------
%get probe spectrum (or pattern, I am not sure yet what NSI gives you)
[E_L2I_AZ_ProbeEx,E_L2I_EL_ProbeEx,E_L2I_AZ_ProbeEy,E_L2I_EL_ProbeEy]...
    = getProbePatterns(TH,PH,0);
% [E_L2I_AZ_ProbeEx,E_L2I_EL_ProbeEx,E_L2I_AZ_ProbeEy,E_L2I_EL_ProbeEy]...
%     = getWR187Patterns(TH,PH,1);

%I rather define the patterns in H,V or X,Y or AZ,EL rather than Co and X
Fx_probeEx = E_L2I_AZ_ProbeEx;
Fy_probeEx = E_L2I_EL_ProbeEx;
Fx_probeEy = E_L2I_AZ_ProbeEy;
Fy_probeEy = E_L2I_EL_ProbeEy;
%Correction - from slater
denominator = Fx_probeEx.*Fy_probeEy - Fx_probeEy.*Fy_probeEx;
F_AZ_corr = (F_AZ.*Fy_probeEy - F_EL.*Fy_probeEx)./denominator;
F_EL_corr = (-F_AZ.*Fx_probeEy + F_EL.*Fx_probeEx)./denominator;

%Spectrum to FF fields (Cos not needed);
EL2IAZ = F_AZ_corr;%.*cosd(TH);
EL2IEL = F_EL_corr;%.*cosd(TH);
warning('No probe correction being used right now.');

EL2IAZ = Fco_Interp;
EL2IEL =  Fcross_Interp;

%Tranform to other polarization definitions
[ETH,EPH] = L2I2ThPh(TH,PH,EL2IAZ,EL2IEL);
[EL3H,EL3V] = ThPh2L3(TH,PH,ETH,EPH);
%  PLOT NF DATA

NFplotter(gridx,gridy,gridco,gridcx);
% PLOT COMPUTED PATTERNS

if pol == 0 %Vertical pol
    EL3co = EL3V;
    EL3x = EL3H;
    %L2I
    EL2Ico = EL2IEL;
    EL2Ix = EL2IAZ;
else
    EL3co = EL3H;
    EL3x = EL3V;    
    %L2I
    EL2Ico = EL2IAZ;
    EL2Ix = EL2IEL;
end

%     [figHandler] = plotPatterns(AZ,EL,EL3co,EL3x,0);
%     suptitle('L3');
%    
%     [figHandler] = plotPatterns(AZ,EL,EL2Ico,EL2Ix,0);
%     suptitle('L2I');
    
    FFexport.AZ = AZ;
    FFexport.EL = EL;
    FFexport.coL3 = EL3co;
    FFexport.crossL3 = EL3x;
    FFexport.coL2I = EL2Ico;
    FFexport.crossL2I = EL2Ix;
%     FFexport.coordinate = 'L3'; % Desired Coordinate System

function [data] = importFFfromfile(filename,mode)
%% Get Co and X from a NSI data file
%mode =1-AzEl(default)/2 - thph/3- Kx/Ky
%data = X,Y,AmpPol1,PhasePol1,AmpPol2,PhasePol2
% For more information, see the TEXTSCAN documentation.
% Ex. filename = '/Users/rodrigolebron/Dropbox/Measurements/NCAR/ARRC Chapter/APAR 8x8/Patterns/1 - Calibrated Uniform/NFCalUnif10.txt';
%     if mode == 1
%         line2find = 'Azimuth (deg)  Elevation (deg)    Amp      Phase   ';
%     elseif mode == 2
%         line2find = 'Theta (deg)  Phi (deg)    Amp      Phase   ';
%     elseif mode == 3
%         line2find = 'Kx (cycles/L)  Ky (cycles/L)    Amp      Phase   ';
%     end
formatSpec = '%f %f %f %f';

%% Open the text file.
fileID = fopen(filename,'r');
 while ~feof(fileID)         
     
    tline = fgetl(fileID);
    
    %look for the first pol
    line2find = 'Probe-1: Lin-';
    IndexPol = strfind(tline,line2find);
    if ~isempty(IndexPol)
        pol = tline(IndexPol+length(line2find))=='0'; %0-V(Ey)/%1-H(Ex)
    end
    
    %look for data
    if mode == 1
        line2find = 'Azimuth (deg)  Elevation (deg)    Amp      Phase   ';
    elseif mode == 2
        line2find = 'Theta (deg)  Phi (deg)    Amp      Phase   ';
    elseif mode == 3
        line2find = 'Kx (cycles/L)  Ky (cycles/L)    Amp      Phase   ';
    end
    IndexC = strfind(tline,line2find);
    if ~isempty(IndexC)
        currentPol = tline(IndexC+length(line2find));
        if currentPol == 'P'
            SData = textscan(fileID,formatSpec,'HeaderLines',0);
            SData = cell2mat(SData);
            [nbr_row, nbr_col] = size(SData);
            nbr_samples = sqrt(nbr_row); %samples along x and y
            %arrange in a 3D matrix
            for i=1:nbr_col
                data(:,:,i) = transp(reshape(SData(:,i),nbr_samples,nbr_samples));
            end
        end
        if currentPol == 'X'
            SData2 = textscan(fileID,formatSpec,'HeaderLines',0);
            Xpolpart = cell2mat(SData2(3:4)); %just get the important part
            [nbr_row,nbr_col] = size(Xpolpart);
            nbr_samples = sqrt(nbr_row); %samples along x and y
            %arrange in a 3D matrix
            for i=5:6
                data(:,:,i) = transp(reshape(Xpolpart(:,i-4),nbr_samples,nbr_samples));
            end
       end
    end
 end
 

 
%rawdata = importdata(filename,' ',49);
fclose(fileID);
% clearvars -except 'data'
end


function [E_L2I_AZ_ProbeEx,E_L2I_EL_ProbeEx,E_L2I_AZ_ProbeEy,E_L2I_EL_ProbeEy] = getProbePatterns(TH,PH,doweplot)

%% Import FF probeEx data in Ludwig2 AZ/EL definition
[Uprobe,Vprobe,Wprobe] = ThPh2UVW(TH,PH);

linewidth = 2;
% filename = '/Users/rodrigolebron/Dropbox/Measurements/NCAR/ARRC Chapter/APAR 8x8/Patterns/1 - Calibrated Uniform/NF/ProbeEx/FFL2AzELKxKyPCProbeEx.txt';%Hpol CalUnif67 data//FFL2ElAzNOPCCalUnif67.txt
%Probe Definition
filename = 'FFL2AzELKxKyPCProbeEx.txt';
dataFF = importFFfromfile(filename,3);
U_meas=dataFF(:,:,1);
V_meas=dataFF(:,:,2);

E_co_DB_file = dataFF(:,:,3);
E_co_Deg_file = dataFF(:,:,4);
E_co_ProbeEx_file = db2mag(E_co_DB_file).*exp(1i*E_co_Deg_file*pi/180);

E_x_DB_file = dataFF(:,:,5);
E_x_Deg_file = dataFF(:,:,6);
E_x_ProbeEx_file = db2mag(E_x_DB_file).*exp(1i*E_x_Deg_file*pi/180);

%PATTERN DATA ----------------
E_L2I_AZ_ProbeEx = interp2(U_meas, V_meas, E_co_ProbeEx_file,Uprobe,Vprobe,'spline');
E_L2I_EL_ProbeEx = interp2(U_meas, V_meas, E_x_ProbeEx_file,Uprobe,Vprobe,'spline');
% TRANSFORM TO OTHER COORDINATES
%Theta & phi
[E_TH_probeEx,E_PH_probeEx]=L2I2ThPh(TH,PH, E_L2I_AZ_ProbeEx, E_L2I_EL_ProbeEx);
%L2AZEL
[E_L3_H_ProbeEx,E_L3_V_ProbeEx]=ThPh2L3(TH,PH, E_TH_probeEx,E_PH_probeEx);
[Ex_probeEx,Ey_probeEx,Ez_probeEx] =ThPh2XYZ(TH,PH, E_TH_probeEx, E_PH_probeEx);

if doweplot
    %for plotting --------------
    F_meas_co = E_L2I_AZ_ProbeEx;
    F_meas_cross = E_L2I_EL_ProbeEx;
    F_comp_co = E_L2I_EL_ProbeEx;
    F_comp_cross = E_L2I_EL_ProbeEx;
    X2plot_meas = Uprobe;
    Y2plot_meas = Vprobe;
    X2plot_comp = Uprobe;
    Y2plot_comp = Vprobe;
    %plot parameters
    ymin = -60;
    ymax = 0;
    Xcut = 0; %where to cut the data
    Ycut = 0; %where to cut the data
    
    %norms
    norm_meas = max(max(F_meas_co));
    norm_Interp = max(max(F_comp_co));
    %Cut in Kx/K or U
    figure();
    subplot(2,1,1);
    hold on;
    indxs = Y2plot_meas==Ycut;
    plot(X2plot_meas(indxs),db(F_meas_co(indxs)/norm_meas),'k');
    indxs = Y2plot_comp==Ycut;
    plot(X2plot_comp(indxs),db(F_comp_co(indxs)/norm_Interp),'g','linewidth',linewidth);
    plot(X2plot_meas(indxs),db(F_meas_cross(indxs)/norm_meas),'k--');
    plot(X2plot_comp(indxs),db(F_comp_cross(indxs)/norm_Interp),'g--','linewidth',linewidth);
    ylim([ymin ymax]);
    xlabel('U');
    ylabel('Magnitude(dB)');
    legend('L2AZ/EL meas','L2AZ/EL from NF');
    grid on;
    grid minor;
    %Cut in Ky/K or V
    subplot(2,1,2);
    hold on;
    indxs = X2plot_meas==Xcut;
    plot(Y2plot_meas(indxs),db(F_meas_co(indxs)/norm_meas),'k');
    plot(Y2plot_meas(indxs),db(F_meas_cross(indxs)/norm_meas),'k--');
    indxs = X2plot_comp==Xcut;
    plot(Y2plot_comp(indxs),db(F_comp_co(indxs)/norm_Interp),'g','linewidth',linewidth);
    plot(Y2plot_comp(indxs),db(F_comp_cross(indxs)/norm_Interp),'g--','linewidth',linewidth);
    ylim([ymin ymax]);
    xlabel('V');
    ylabel('Magnitude(dB)');
    grid on;
    grid minor;
    suptitle('L3 FF - Probe Ex');
    
    
    %plot cartesian components
%    plotPatterns(AZ,EL,E_L2I_AZ_ProbeEx,E_L2I_EL_ProbeEx,[0])
%    suptitle('Probe Ex- L2 AZ/EL')
    % figure
    % surf(U,V,db(Ex_probeEx))
    % shading interp
    % view(2)
    % xlabel('U')
    % ylabel('V')
    % title('Ex - Probe Ex')
    % %Ey
    % figure
    % surf(U,V,db(Ey_probeEx))
    % shading interp
    % view(2)
    % xlabel('U')
    % ylabel('V')
    % title('Ey - Probe Ex')
end
%% Import FF probe data Ey
% filename = '/Users/rodrigolebron/Dropbox/Measurements/NCAR/ARRC Chapter/APAR 8x8/Patterns/1 - Calibrated Uniform/NF/ProbeEy/FFL2AzELKxKyPCProbeEy.txt';%Hpol CalUnif67 data//FFL2ElAzNOPCCalUnif67.txt
filename = 'FFL2AzELKxKyPCProbeEy.txt';
dataFF = importFFfromfile(filename,3);
U_meas=dataFF(:,:,1);
V_meas=dataFF(:,:,2);

E_co_DB_file = dataFF(:,:,3);
E_co_Deg_file = dataFF(:,:,4);
E_co_ProbeEy_file = db2mag(E_co_DB_file).*exp(1i*E_co_Deg_file*pi/180);

E_x_DB_file = dataFF(:,:,5);
E_x_Deg_file = dataFF(:,:,6);
E_x_ProbeEy_file = db2mag(E_x_DB_file).*exp(1i*E_x_Deg_file*pi/180);

%PATTERN DATA
E_L2I_AZ_ProbeEy = interp2(U_meas, V_meas, E_x_ProbeEy_file,Uprobe,Vprobe,'spline');
E_L2I_EL_ProbeEy = interp2(U_meas, V_meas, E_co_ProbeEy_file,Uprobe,Vprobe,'spline');
% TRANSFORM TO OTHER COORDINATES
%Theta & phi
[E_TH_probeEy,E_PH_probeEy]=L2I2ThPh(TH,PH, E_L2I_AZ_ProbeEy, E_L2I_EL_ProbeEy);
%L2AZEL
[E_L3_H_ProbeEy,E_L3_V_ProbeEy]=ThPh2L3(TH,PH, E_TH_probeEy,E_PH_probeEy);
[Ex_probeEy,Ey_probeEy,Ez_probeEy] =ThPh2XYZ(TH,PH, E_TH_probeEy, E_PH_probeEy);
%THE ACTUAL CORRECTION FACTORS

if doweplot
    %for plotting --------------
    F_meas_co = E_L2I_EL_ProbeEy;
    F_meas_cross = E_L2I_AZ_ProbeEy;
    F_comp_co = E_L2I_AZ_ProbeEy;
    F_comp_cross = E_L2I_AZ_ProbeEy;
    X2plot_meas = Uprobe;
    Y2plot_meas = Vprobe;
    X2plot_comp = Uprobe;
    Y2plot_comp = Vprobe;
    %plot parameters
    ymin = -60;
    ymax = 0;
    Xcut = 0; %where to cut the data
    Ycut = 0; %where to cut the data
    
    %norms
    norm_meas = max(max(F_meas_co));
    norm_Interp = max(max(F_comp_co));
    
    %Cut in Kx/K or U
    figure();
    subplot(2,1,1);
    hold on
    indxs = Y2plot_meas==Ycut;
    plot(X2plot_meas(indxs),db(F_meas_co(indxs)/norm_meas),'k');
    indxs = Y2plot_comp==Ycut;
    plot(X2plot_comp(indxs),db(F_comp_co(indxs)/norm_Interp),'g','linewidth',linewidth);
    plot(X2plot_meas(indxs),db(F_meas_cross(indxs)/norm_meas),'k--');
    plot(X2plot_comp(indxs),db(F_comp_cross(indxs)/norm_Interp),'g--','linewidth',linewidth);
    ylim([ymin ymax]);
    xlabel('U');
    ylabel('Magnitude(dB)');
    legend('L2AZ/EL meas','L2AZ/EL from NF');
    grid on;
    grid minor;
    %Cut in Ky/K or V
    subplot(2,1,2);
    hold on;
    indxs = X2plot_meas==Xcut;
    plot(Y2plot_meas(indxs),db(F_meas_co(indxs)/norm_meas),'k');
    plot(Y2plot_meas(indxs),db(F_meas_cross(indxs)/norm_meas),'k--');
    indxs = X2plot_comp==Xcut;
    plot(Y2plot_comp(indxs),db(F_comp_co(indxs)/norm_Interp),'g','linewidth',linewidth);
    plot(Y2plot_comp(indxs),db(F_comp_cross(indxs)/norm_Interp),'g--','linewidth',linewidth);
    ylim([ymin ymax]);
    xlabel('V');
    ylabel('Magnitude(dB)');
    grid on;
    grid minor;
    suptitle('L3 FF - Probe Ey ');

    %plot cartesian components
    %Ex
    %plot cartesian components
%    plotPatterns(AZ,EL,E_L2I_AZ_ProbeEy,E_L2I_EL_ProbeEy,[0])
%    suptitle('Probe Ey - L2 AZ/EL')
    % figure
    % surf(U,V,db(Ex_probeEy))
    % shading interp
    % view(2)
    % xlabel('U')
    % ylabel('V')
    % title('Ex - Probe Ey')
    % %Ey
    % figure
    % surf(U,V,db(Ey_probeEy))
    % shading interp
    % view(2)
    % xlabel('U')
    % ylabel('V')
    % title('Ey - Probe Ey')
end
end
function [figHandler] = plotPatterns(AZ,EL,PatternCmplx_Co,PatternCmplx_X,cutAngleOptions,autoscale)
%% PLOTTING

if nargin >= 5
    options = cutAngleOptions;
else
    options = -45:5:45;
end

if nargin <= 5
    autoscale = 1;

end


Pattern2plot_Co = PatternCmplx_Co;
Pattern2plot_X = PatternCmplx_X;

%get norm (MAX) val
normCo = max(max(Pattern2plot_Co));
normX = max(max(Pattern2plot_X));
norm = max([normCo normX]);
if normCo > normX
    index2norm = Pattern2plot_Co==norm;
else
    index2norm = Pattern2plot_X==norm;
end

%2D plot - Co
figure();
subplot(2,2,3);
surf(AZ,EL,db(Pattern2plot_Co/norm));
colorbar
shading interp;
axis square;
caxis([-50 0]);
if autoscale
    xlim([-65 65]);
    ylim([-65 65]);
end
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
view(2);

%2D plot - X
subplot(2,2,2);
surf(AZ,EL,db(Pattern2plot_X/norm));
% title('L3 FF Cross-Polarization');
colorbar;
shading interp;
axis square;
caxis([-50 0]);
if autoscale
    xlim([-65 65]);
    ylim([-65 65]);
end
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
view(2);




%fetch for cut angles

normElAngle = EL(index2norm);
normAzAngle = AZ(index2norm);
[angleErrorEl, index2ElAngle] = min(abs(options-normElAngle));
[angleErrorAz, index2AzAngle] = min(abs(options-normAzAngle));
targetElAngle = options(index2ElAngle);
targetAzAngle = options(index2AzAngle);

%AZ cut
[~, index2EL]= min(abs(EL(:)-targetElAngle));
index2plot  = EL == EL(index2EL(1));
h1 =subplot(2,2,1);
hold on
plot(AZ(index2plot),db(Pattern2plot_Co(index2plot)/norm));
plot(AZ(index2plot),db(Pattern2plot_X(index2plot)/norm),'--');
if autoscale
    xlim([-65 65]);
end
ylim([-50 0]);
xlabel('Azimuth (deg)');
ylabel('Amplitude(dB)');
grid on;
grid minor;
axis square;
%ylim([-65 65])

%EL cut
[~, index2Az]= min(abs(AZ(:)-targetAzAngle));
index2plot  = AZ == AZ(index2Az(1));
h1 =subplot(2,2,4);
%set(h1, 'Ydir', 'reverse')
%set(h1, 'YAxisLocation', 'Right')
hold on
plot(db(Pattern2plot_Co(index2plot)/norm),EL(index2plot));
plot(db(Pattern2plot_X(index2plot)/norm),EL(index2plot),'--');
if autoscale
    ylim([-65 65]);
end
xlim([-50 0]);
ylabel('Elevation (deg)');
xlabel('Amplitude(dB)');
grid on;
grid minor;
axis square;

figHandler = gcf;
end
function [U,V,W] = AzEl2UVW(AZ,EL)
U = cosd(EL).*sind(AZ);
V = sind(EL);
W = cosd(EL).*cosd(AZ);
end
function [TH,PH] = UVW2ThPh(U,V,W)
PH = atan2(V,U)*180/pi;
TH = asind(V./sind(PH));
%issue when PH = 0, then 
index = isnan(TH);%
TH(index) = acosd(W(index));
%TH(index) = abs(AZ(index));
end
function [U,V,W] = ThPh2UVW(TH,PH)
U = sind(TH).*cosd(PH);
V = sind(TH).*sind(PH);
W = cosd(TH);
end
function [ETH,EPH] = L2I2ThPh(TH,PH,EL2IAZ,EL2IEL)
norm = 1./sqrt(1-(sind(TH)).^2.*(sind(PH)).^2);
ETH = (cosd(PH).*EL2IAZ         +cosd(TH).*sind(PH).*EL2IEL).*norm;
EPH  = (-cosd(TH).*sind(PH).*EL2IAZ         +cosd(PH).*EL2IEL).*norm;
end
function [EL3co,EL3x] = ThPh2L3(TH,PH,ETH,EPH)
EL3co = cosd(PH).*ETH -sind(PH).*EPH;
EL3x  = sind(PH).*ETH +cosd(PH).*EPH;
end
function [EX,EY,EZ] = ThPh2XYZ(TH,PH,ETH,EPH)
EX =  cosd(TH).*cosd(PH).*ETH   -sind(PH).*EPH;
EY =  cosd(TH).*sind(PH).*ETH   +cosd(PH).*EPH;
EZ = -sind(TH).*ETH;
end
function [] = NFplotter(gridx,gridy,gridco,gridcx)
%Feel free to remove the cross pol when the time comes.
%this function is fairly easy, just send in the x and y grids along with
%the data to be plotted. It will make a 3d surf plot for you.
    figure();
    surf(gridx,gridy,gridco- max(max(gridco)));
    shading interp
    colorbar
    view(2)
    xlabel('X(m)')
    ylabel('Y(m)')
    title('NF Magnitude - Co')
    h = colorbar();
    ylabel(h,'Magnitude (dB)');
%     saveas(gcf,'CoFF.png');

    figure();
    surf(gridx,gridy,gridcx- max(max(gridco)));
    shading interp
    colorbar
    view(2)
    xlabel('X(m)')
    ylabel('Y(m)')
    title('NF Magnitude - Cross')
    h = colorbar();
        ylabel(h,'Magnitude (dB)');
    
    figure();
    surf(gridx,gridy,gridpco);
    shading interp
    colorbar
    view(2)
    xlabel('X(m)')
    ylabel('Y(m)')
    title('NF Co Phase')
    h = colorbar();
        ylabel(h,'Phase (deg)');

    figure();
    surf(gridx,gridy,gridpcx);
    shading interp
    colorbar
    view(2)
    xlabel('X(m)')
    ylabel('Y(m)')
    title('NF Cross Phase')
    h = colorbar();
        ylabel(h,'Phase (deg)');
%     saveas(gcf,'Phaseco.png');
    
    end
end
