clear all
%% Parameters
del = 19*pi/180; %dip angle in radians
psi = pi/4-del; % sigma 1 plunge angle in radians
gamma = 2.1e-11;% Shear strain rate 66 mm/yr over 100 m thickness
e = gamma/(3^0.5); %Conversion from simple shear to equaivalent axial-symmetric
vel = 2.1e-9; % Shear velocity 66 mm/yr
velref = 1e-6; % Reference velocity for frictional deformation
L = 0.4; %Pore pressure factor, lambda
M = readmatrix('Shumagin_EM.xlsx');
Depth = -1*M(:,4); %Depth in km
TK = M(:,2); %Temperature in K
TC = M(:,1); %Temperature in C
Pressure =  M(:,3)/10; %pressure in MPa (divide by 10 bc in bars)
Pf = L*Pressure; %or = M(:,5);
%% Constants
R = 8.314;
%% Lawsonite Blueshist Friction from Sawai et al., 2016
abT = 0.008-2.3e-4*(TC-200); %a-b as a fuction of temperature between 200 and 300 C
abT(TC<200) = 0.008; %a-b below 200 C
abT(TC>300) = -0.015; %a-b above 300 C
mutT = 0.8-0.003*(TC-200); %friction coefficient as a fuction of temperature between 200 and 300 C
mutT(TC<200) = 0.8; %friction coefficient below 200 C
mutT(TC>300) = 0.5; %friction coefficient above 300 C
mutTr = mutT+abT*log(vel/velref); %friction coefficient corrected for velocity
ttT = (mutTr.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mutTr*(cos(2*del+2*psi)-cos(2*psi))); %Shear stress in MPa
%% Hufford pre preprint 2024
A = 2.23e5;
Q = 341e3; %J/mol
n = 3;
sigma_MPa_dc = (e./(A*exp(-Q./(R*TK)))).^(1/n); %differential stress from flow law in MPa
tau_dc = sigma_MPa_dc/(3^0.5); %differential stress converted to equivalent simple shear stress in MPa
%% Total Deformation
total = (ttT.^(-20)+tau_dc.^(-20)).^(-1/20);
%% Figures
figure(1)
hold on
plot(ttT(TC>200), -1*Depth(TC>200), tau_dc(TC>200), -1*Depth(TC>200), total(TC>200), -1*Depth(TC>200))
box on
xlim([10^(-1) 1000])
ylim([-80 -10])
xlabel 'Shear Stress (MPa)'
ylabel 'Depth (km)'
%
%data = [ttT(TC>200) tau_dc(TC>200)]; %
%writematrix(data,'Shumagin_EM.xlsx','Range','T:U')
