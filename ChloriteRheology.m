clear
%% Parameters
del = 19*pi/180; %dip angle in radians
psi = pi/4-del; % sigma 1 plunge angle in radians
vel = 2.1e-9;% plate velocity 66 mm/yr
velref = 1e-6; %reference velocity for frictional deformation
L = 0.4; %pore pressure factor, lambda
M = readmatrix('Shumagin_EM.xlsx');
Depth = -1*M(:,4); %Depth in km
TK = M(:,2); %Temperature in K
TC = M(:,1); %Temperature in C
Pressure =  M(:,3)/10; %pressure in MPa (divide by 10 bc in bars)
Pf = L*Pressure; %or = M(:,5);
%% Chlorite friction using Okamato et al., 2019
% %
abCM = 0.005; %a-b
mutCM = zeros(size(TC)); %create friction vector
mutCM(TC<400) = 0.275; %friction coefficient less than 400 C
mutCM(TC==400) = 0.275; %friction coefficient equal to 400 C
mutCM(TC>400) = 0.275+7.5e-4*((TC(TC>400)-400)); %friction coefficient 400 to 500 C
mutCM(TC>500) = 0.35; %friction coefficient greater than 500 C
mutCMR = mutCM +abCM*log(gamma./gammaref);
ttCM = (mutCMR.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mutCMR*(cos(2*del+2*psi)-cos(2*psi)));
%% Figures
figure(1)
hold on
plot(ttCM, -1*Depth)
box on
xlim([0 500])
ylim([-80 -10])
xlabel 'Shear Stress (MPa)'
ylabel 'Depth (km)'
%
%writematrix(ttCM,'Shumagin_EM.xlsx','Range','R1')