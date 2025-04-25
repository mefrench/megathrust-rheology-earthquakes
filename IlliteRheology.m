clear all
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
%% constants
R = 8.314; %gas constant, (K*J/mol);
%% Illite friction Den Hartog et al., 2012
abT = 0.053-2e-4*TC; %temperature dependence of a-b
abT(TC<201) = 0.01; %a-b below 200 C
mutT(TC<301) = 0.6+3.3e-4*(TC(TC<301)-150); %friction coefficient below 300 C
mutT(TC>=301) = 0.65+0.003*((TC(TC>301)-300)); %friction coefficient from 300 to 350 C
mutT(TC>350) = 0.8+0.001*((TC(TC>350)-350)); %friction coefficient above 350 C
%
mutTr = mutT+abT*log(vel/velref); %friction coefficient adjusted for rate dependence
ttT = (mutTr.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mutTr*(cos(2*del+2*psi)-cos(2*psi))); %shear stress in MPa
%% Figures
figure(1)
hold on
plot(ttT(TC<300), -1*Depth(TC<300))
xlim([10^(-1) 500])
ylim([-80 -10])
xlabel 'Shear Stress (MPa)'
ylabel 'Depth (km)'
%
%writematrix(ttT(TC<300),'Shumagin_EM.xlsx','Range','N1')