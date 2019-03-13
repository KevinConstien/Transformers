function [] = kplot(FF)
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
if sum(fields) < length(vars);
    error('Incorrect input. See the required inputs for more info.')
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

%Let's made some 3D plots;

%Co
figure();
% subplot(2,2,1);
surf(AZ,EL,coNorm);
title('Co-pol');
setup3D();

%Cross
figure();
% subplot(2,2,3);
surf(AZ,EL,cxNorm);
title('Cross-pol');
setup3D();
[cmin, ~] = caxis;
caxis([cmin, 0]);

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
    warning(mywarning);
    if length(row) + length(col) == 0
        error('No maximum value has been found to place the Az/EL cuts on.');
    end
end

lw = 1;

%Horizontal Plane
figure();
% subplot(2,2,2);
plot(AZ(row,:),coNorm(row,:),'LineWidth',lw); hold on;
plot(AZ(row,:),cxNorm(row,:),'-','LineWidth',lw);
title("Radiation Pattern Co Horizonatal Plane (\phi = 0^{\circ}) " + string(FF.coordinate));
% xlim([min(min(AZ(row,:))),max(max(AZ(row,:)))]);
ylim([-50 0]);
xlim([-40 40]);

setup2D();
xlabel('Azimuth (deg)');

%Vertical Plane
figure();
% subplot(2,2,4);
plot(EL(:,col),coNorm(:,col),'LineWidth',lw);hold on;
plot(EL(:,col),cxNorm(:,col),'-','LineWidth',lw);
title("Radiation Pattern Co Vertical Plane (\phi = 90^{\circ}) " + string(FF.coordinate));
% xlim([min(EL(:,col)),max(EL(:,col))]);
ylim([-50 0]);
xlim([-40 40]);

setup2D();
xlabel('Elevation (deg)');
if isfield(FF,'title')
    mytitle = FF.title + " "+FF.coordinate;
    h = suptitle(mytitle);
    set(h,'FontSize',24,'FontWeight','normal');
    
end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,0,1,.96]);%[0, 0.04, 1, 0.96]);

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
        legend('Co','Cross','Location','Southeast')
        ylabel('Magnitude(dB)');
        grid minor;
        axis normal;        
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


