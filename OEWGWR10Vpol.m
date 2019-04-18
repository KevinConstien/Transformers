%%This code outputs a Ludwig 2I definition, H/V components in Azimuth/Elevation coordinates for 
%an OEWG polarized along the y axis(elevation/vertical)
%Produces approximate results, as WR-10 is significantly smaller than other
%OEWG cited in paper
%[AZ,EL,EAZ,EEL] = OEWGWR10Vpol(freq)  
function [AZ,EL,EAZ,EEL] = OEWGWR10Vpol(freq) 
%set(0,'defaultTextInterpreter','tex'); %trying to set the default


f = freq;
c = 3e8;
lambda = c/f;
k = (2*pi)/lambda;
Z0 = 377;
a = 2.54e-3; %%Dimensions for WR-10 Waveguide
b = 1.27e-3;
E0 = 1;
A_e = 1;
gamma = 0;
theta = (-80:1:80).*((2*pi)/(360));
phi = (0:1:179).*((2*pi)/(360));
[TH,PH] = meshgrid(theta,phi);
%transform to AZ/EL
[U,V,W] = ThPh2UVW(TH*180/pi,PH*180/pi);
[AZ,EL] = UVW2AzEl(U,V,W);
%phi = 45;
r = 1000;

%A_h = (-1*i*(k^2)*a*b*E0)/8;

Bk = sqrt(1-(pi/(k*a))^2); %Beta/k calculation

x = ((k*b)/2).*sin(TH); %Calculation variable

Ee_num = A_e.*(1 + Bk.*cos(TH)+gamma*(1 - Bk.*cos(TH))).*sin(x); 
Ee_den = (1+Bk+gamma*(1 - Bk)).*x;

Ee = Ee_num./Ee_den;
EedB = db(Ee);
EeNormdB = EedB - max(max(EedB));

y = ((k*a)/2).*sin(TH);

Eh = A_e.*((pi/2)^2).*cos(TH).*(cos(y)./(((pi/2)^2)-(y.^2)));
EhdB = db(Eh);
EhNormdB = EhdB - max(max(EhdB));

z = (exp(1i*k*r))/(k*r);

Eth = Ee.*sin(PH); Eph = Eh .*cos(PH);
[EAZ,EEL] = ThPh2L3(TH*180/pi,PH*180/pi,Eth,Eph);
E = db(z.*(Ee.*sin(PH) + Eh.*cos(PH)));
Enorm = E - max(E);
H = db((1/Z0).*z.*(Ee.*sin(PH) - Eh.*cos(PH)));
Hnorm = H - max(H);


th_plot = theta*(360)./(2*pi);
%plot(th_plot,Ee,'linewidth',3); hold on;
%ylim([-60 0]);
% grid on;
% grid minor;
% ylabel('Amplitude (dB)');
% xlabel('Theta (deg)');

% %plot in 3d - Eth & Eph
% figure;surf(AZ,EL,db(Eph));%surf(TH,PH,db(Ephi));
% shading interp
% view(2)
% xlabel('Azimuth')%'Theta (rad)')
% ylabel('Elevation')%'Phi (rad)')
% title('E_{ph}')
% %plot in 3d
% figure;surf(AZ,EL,db(Eth));%surf(TH,PH,db(Ephi));
% shading interp
% view(2)
% xlabel('Azimuth')%'Theta (rad)')
% ylabel('Elevation')%'Phi (rad)')
% title('E_{th}')
% 
% %plot in 3d - Eaz & Eth
% figure;surf(AZ,EL,db(EAZ));%surf(TH,PH,db(Ephi));
% shading interp
% view(2)
% xlabel('Azimuth')%'Theta (rad)')
% ylabel('Elevation')%'Phi (rad)')
% title('E_{az}')
% %plot in 3d
% figure;surf(AZ,EL,db(EEL));%surf(TH,PH,db(Ephi));
% shading interp
% view(2)
% xlabel('Azimuth')%'Theta (rad)')
% ylabel('Elevation')%'Phi (rad)')
% title('E_{El}')
% 
% %% cut EV and EH
% 
% figure;
% subplot(2,1,1)
% plot(AZ(EL==0),db(EEL(EL==0)))
% ylim([-10 0])
% xlabel('Azimuth')
% grid on 
% grid minor
% subplot(2,1,2)
% plot(EL(AZ==0),db(EEL(AZ==0)))
% ylim([-10 0])
% xlabel('Elevation')
% grid on
% grid minor

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
%[AZ,EL] = UVW2AzEl(U,V,W)
%U,V,W  meshgrid type
%AZ EL defined by Masters et al. 
function [AZ,EL] = UVW2AzEl(U,V,W)
EL = asind(V);
AZ = asind(U./cosd(EL));
end
end