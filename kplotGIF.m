function [] = kplotGIF(FF)
% KPLOT A far-field data plotter.
%   kplot(FF) plots the far-field data specified in the FF struct. Make
%   sure the FF input contains the fields FF.AZ, FF.EL, FF.crossL3, FF.coL2I,
%   and FF.crossL2i. If one of these fields is missing, the program will
%   not run. If using this function along with Kevin Constien's NF2FF
%   function, then the appropriate struct will be output. To select the
%   coordinate system set FF.coordinate to either 'L3' or 'L2.' The default
%   is L3. 
%   
%   Use FF.title = ('My Title') to set a title to the whole figure.
%
%   NOTE: Each of these variables are 2D grids. AZ and EL are 2D so that a 
%   See also NF2FFFUNCTION, EXAMPLESCRIPT.


% clc
% close all
% clear all
% load('FF.mat');

% Checking to make sure correct input.
vars = {'AZ','EL','coL3','crossL3','coL2I','crossL2I'};
fields = isfield(FF,vars);
if sum(fields) < length(vars)
    error('Incorrect input. See the required inputs for more info.')
end

vars = {'azangle','elangle'};
if sum(isfield(FF,vars)) == 2   
    elangle = FF.elangle;
    azangle = FF.azangle;
else
    elangle = 0;
    azangle = 0;
end


nocoordinate = "No coordinate system has been chosen. Defaulting to Ludwig 3.";
invalidcoor = "FF.coordinate is not set to a valid coordinate system. " + ...
    "Defaulting to Ludwig 3.";
if ~isfield(FF,'coordinate');
    warning(nocoordinate);
    FF.coordinate = 'L3';
    co = FF.coL3;
    cross = FF.crossL3;
else
    
    switch FF.coordinate
        case {'L3'}
            co = FF.coL3;
            cross = FF.crossL3;
            
        case{'L2'}
            co = FF.coL2I;
            cross = FF.crossL2I;
        otherwise
            warning(invalidcoor);
            co = FF.coL3;
            cross = FF.crossL3;
    end
end

AZ = FF.AZ;
EL = FF.EL;


%Normalizeing to 0 dB.
maxCo = max(max(db(co)));
maxCx = max(max(db(cross)));
globalMax = max(maxCo,maxCx);
coNorm = db(co) - globalMax;
cxNorm = db(cross) - globalMax;

% coNorm = db(co);
% cxNorm = db(cross);
set(0,'defaultAxesFontSize',17);
set(0,'defaultAxesFontName','Times');

Hfig = .8;
Wfig = .8;
figHandler = figure('Units','normalized','Position',[0.2 0.1 Wfig Hfig]);
set(0, 'currentfigure', gcf);
set(figHandler,'Units','pixel')
figPos = get(figHandler,'Position');
set(figHandler,'Position', ...
    [figPos(1) figPos(2) figPos(4) figPos(4)])

 %DRAW AXES
H1 = 0.45;
W1 = 0.45;
X1 = (1-W1)/2+0.2;
Y1 = (1-H1)/2-0.2;
ax1 = axes('Parent',figHandler,...
    'Units','normalized','position',[X1 Y1 W1 H1]);
% set(ax1,'Units','pixel')
pos1 = get(ax1,'Outerposition');


%AZ cut - to the top
H2 = (1 - H1)/2;
W2 = W1;
X2 = X1;
Y2 =H1+Y1+0.08;
ax2 = axes('Parent',figHandler,...
    'Units','normalized','position',[X2 Y2 W2 H2]);
box on

%EL cut - to the left
H3 = H1;
W3 =(1 - W1)/2;
X3 = X1-W3-.08;
Y3 = Y1;
ax3 = axes('Parent',figHandler,...
    'Units','normalized','position',[X3 Y3 W3 H3]);
set(ax3, 'Xdir', 'reverse')
ax3.YTickLabelRotation = 90;
set(ax3, 'YAxisLocation', 'Right')
box on

%top left corner
% H4 = H2;
% W4 = W3;
% X4 = X3-0.05;
% Y4 = Y2+0.05;
% ax4 = axes('Parent',figHandler,...
%     'Units','normalized','position',[X4 Y4 W4 H4]);

%%%%%%%%%%%%%%%%%%%%%%%
%Let's made some 3D plots;
%3D 3D 3D
%%%%%%%%%%%%%%%%%%%%%%%
%Co
% figure();
% subplot(2,2,4);
axes(ax1)
surf(AZ,EL,coNorm);
% title('Co-pol');
setup3D();
xlim([-50 50]);
ylim([-50 50]);
caxis([-70 0]);
set(gcf,'color','w');
% %Cross
% figure();
% % subplot(2,2,3);
% surf(AZ,EL,cxNorm);
% title('Cross-pol');
% setup3D();
% [cmin, ~] = caxis;
% caxis([cmin, 0]);
% 
mywarning = "The cross-pol data has a higher maximum! "...
    + "The co and cross data may be switched.";
[row,col] = find(db(co) == globalMax);

%This find the highest point in the data and sets it as the center. It is
%then used to select the principle plane in AZ and El. This can be changed
%easily. Simply select the row and column that you would like to make the
%desired cut. I have not added the functionality to plot the D-plane yet,
%but if you want to add that, then have fun. 

if length(row) + length(col) == 0
    [row,col] = find(db(cross) == globalMax);
%     warning(mywarning);
    if length(row) + length(col) == 0
        error('No maximum value has been found to place the Az/EL cuts on.');
    end
end

lw = 2;
%%%%%%%%%%%%%%%%
% 2D cuts
%%%%%%%%%%%%%%%%
%Horizontal Plane
axes(ax2);
% [row,~] = (find(FF.EL == elangle)); row = unique(row);
plot(AZ(row,:),coNorm(row,:),'LineWidth',lw); hold on;
% plot(AZ(row,:),cxNorm(row,:),'-','LineWidth',lw);
% title("Radiation Pattern Co Horizonatal Plane (EL angle = "+ elangle +...
%     "^{\circ}) ");
setup2D();
xlim([-50, 50]);
ylim([-40,0]);
ylabel('Magnitude (dB)');



%Vertical Plane
% mag77 = double(flipud(["7.582718";"5.7151633";"5.1994396";"9.521879";"15.618265";"19.844564";"21.979406"...
%     ;"24.139982";"27.675507";"30.938269";"33.901056";"36.827619";"38.535035";"42.032949";...
%     "42.976688";"38.045688";"35.261756";"34.284386";"33.782255";"33.941335";"32.920417";...
%     "32.216306";"27.433687";"23.561932";"22.161476";"20.082277";"19.037522";"17.773773";...
%     "15.914499";"13.909487";"12.223983"]));
% 
% mag77 = mag77 - max(mag77);
% 
% elevations = double(flipud(["15";"14";"13";"12";"11";"10";"9";"8";"7";"6";"5";"4";"3";"2";...
%     "1";"-0";"-1";"-2";"-3";"-4";"-5";"-6";"-7";"-8";"-9";"-10";"-11";"-12";...
%     "-13";"-14";"-15"]));
axes(ax3);
% [~,col] = (find(FF.AZ == azangle)); col = unique(col);
% close gcf;
% figure();
% plot(EL(:,col),coNorm(:,col),'LineWidth',lw);hold on;
% plot(elevations,mag77,'LineWidth',lw);
plot(coNorm(:,col),EL(:,col),'LineWidth',lw);hold on;
% plot(mag77,elevations,'LineWidth',lw);
% legend('Measured','Metawave Data');
% plot(EL(:,col),cxNorm(:,col),'-','LineWidth',lw);
% title("Radiation Pattern Co Vertical Plane (AZ angle = "+ azangle +...
%     "^{\circ})");
setup2D();
xlim([-40,0]);
ylim([-50,50]);
set(ax3, 'Xdir', 'reverse')
ax3.YTickLabelRotation = 90;
set(ax3, 'YAxisLocation', 'Right')
xlabel('Magnitude (dB)');
yticks([-50,50]);


set(0,'defaultAxesFontSize',11);
set(0,'defaultAxesFontName','Helvetica');
% ylabel('Elevation Angle (deg)');
% title('Elevation Cut 77GHz');

% xlabel('Elevation (deg)');
% if isfield(FF,'title')
%     mytitle = FF.title + " "+FF.coordinate;
%     h = suptitle(mytitle);
%     set(h,'FontSize',24,'FontWeight','normal');
%     
% end

% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,0,1,.96]);%[0, 0.04, 1, 0.96]);

% %Polar Plot... wouldnt recommend using this.
% figure();
% subplot(1,2,1);
% polarplot(deg2rad(AZ(row,:)) + (pi/2),coNorm(row,:),'LineWidth',2);
% rlim([-60 0]);
% title('Horizontal Plane');
% thetalim([0,180]);
% 
% subplot(1,2,2);
% polarplot(deg2rad(EL(:,col)) + (pi/2),coNorm(:,col),'LineWidth',2);
% rlim([-60 0]);
% title('Vertical Plane');
% thetalim([0,180]);


% figure();
% plot(deg2rad(AZ(row,:)),coNorm(row,:));
% rlim = [-40 0];
% xlabel('Azimuth (deg)');
% ylabel('Magnitude (dB)');

    function [] = setup2D()
%         legend('Co','Cross','Location','Southeast')
%         ylabel('Magnitude(dB)');
        grid minor;
        axis normal;
%         ylim([-35 0]);
    end


    function [] = setup3D()
        xlabel('Azimuth (deg)');
        ylabel('Elevation (deg)');
        zlabel('Magnitude (dB)');
        axis normal;
        h = colorbar();
        ylabel(h,'Magnitude (dB)');
        colormap jet;
        shading interp;
        grid off;
        view(2);
    end
end


