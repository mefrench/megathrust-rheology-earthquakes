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
%% Talc friction Moore and Lockner, 2008
abT = -1.81e-5*TC+0.00915; % a-b
mutT = 0.1744-0.000157*TC; % temperature dependent friction coefficient
mutTr = mutT+abT*log(vel/velref); % adjust friction coefficent for rate dependence 
ttT = (mutTr.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mutTr*(cos(2*del+2*psi)-cos(2*psi))); %shear stress in MPa
%% Semi-brittle Talc Escartin et al., 2008
muint_e = 0.139-0.0001633*(TC-23); %internal friction as a function of temperature
tau_0 = 13; % cohesion at room T and 400 C
taue_ei = (muint_e.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+muint_e*(cos(2*del+2*psi)-cos(2*psi))); %shear stress in MPa from internal friction
tau_ei = tau_0+taue_ei; %total shear stress in MPa
%% Boneh et al., (2023)
muB0 = 0.13; %friction coefficient up to 400 C
muBT = muB0-4e-4*(TC-400); %Temperature dependence from 400 to 700C
muBtlc = zeros(size(TC)); %create friction vector
muBtlc(TC<400) = muB0; %friction coefficient up to 400 C
muBtlc(TC==400) = muB0; %friction coefficient up to 400 C
muBtlc(TC>400) = muBT(TC>400); %friction coefficient from 400 to 700 C
tB = (muBtlc.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+muBtlc*(cos(2*del+2*psi)-cos(2*psi)));
%% Total Talc Deformation
tauT = (ttT.^(-20)+tau_ei.^(-20)).^(-1/20); %Moore and Lockner and Escartin
%% Figures
figure(1)
hold on
plot(tauT, -1*Depth, ttT, -1*Depth, tau_ei, -1*Depth)
box on
xlim([0 500])
ylim([-80 -10])
xlabel 'Shear Stress (MPa)'
ylabel 'Depth (km)'
%
%writematrix((ttT.^(-20)+tau_ei.^(-20)).^(-1/20),'Shumagin_EM.xlsx','Range','S1')
