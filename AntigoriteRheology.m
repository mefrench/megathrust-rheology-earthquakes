%% Parameters
del = 19*pi/180; %dip angle in radians
psi = pi/4-del; % sigma 1 plunge angle in radians
gamma = 4.2e-12;% Shear strain rate 66 mm/yr over 500 m thickness
e = gamma/(3^0.5); %Conversion from simple shear to equaivalent axial-symmetric
vel = 2.1e-9; % Shear velocity (m/s) 66 mm/yr
velref = 1e-6; % Reference velocity (m/s) for frictional deformation
L = 0.4; %Pore pressure factor, lambda
M = readmatrix('Shumagin_EM.xlsx');
Depth = -1*M(:,4); %Depth in km
TK = M(:,2); %Temperature in K
TC = M(:,1); %Temperature in C
Pressure =  M(:,3)/10; %pressure in MPa (divide by 10 bc in bars)
Pf = L*Pressure; %or = M(:,5);
%% Constants
R = 8.314;
%% Lizardite from Moore and Lockner, 1997
abl(TC<200) = -0.003; %a-b for less than 200 C
abl(TC>200) = 0.005; %a-b for greater than 200 C
mulT(TC<200) = 0.52+8e-4*(TC(TC<200)-200); %temperature dependent friction for less than 200 C
mulT(TC>=200) = 0.52; %friction at greater than 200 C
mul = mulT+abl*log(vel/velref); %friction adjusted for rate dependence
ttl = (mul.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mul*(cos(2*del+2*psi)-cos(2*psi))); %shear stress (MPa)
%ttl(TC>400)=1000;
%% Antigorite friction Takahashi et al., 2011
Tab = [25 100 200 300 350 400 450]; %temperatures at which a-b was measured
abT = [0.00915 0.01078 0.00818 0.00629 0.01609 0.01835 -0.00227]; %measurements of a-b
abs = interp1(Tab', abT', TC, 'linear',-0.00227); %linear interpolation between measurements of a-b
muT = [0.55 0.53 0.52 0.52 0.52 0.65 0.70 0.70]; % measurements of friction coefficient
Tmu = [23 100 200 300 400 450 500 600]; % temperatures at which measurements of friction were made
musi = interp1(Tmu', muT', TC, 'linear','extrap'); %linear interpolation of friction between temperatures
mus = musi +abs*log(vel/velref); %friction coefficient adjusted for rate-dependence
tts = (mus.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mus*(cos(2*del+2*psi)-cos(2*psi))); %shear stress (MPa)
%% Antigorite Rheology from Proctor and Hirth (2016)
abb = 0.01; %a-b
mubT = -0.0008*TC+0.4633; %internal friction as a function of temperature up to 500 C
mubT(TC>500) = 0.07; %internal friction above 500 C
mub = mubT+ abb*log(vel/velref); %internal friction adjusted for rate
tfp = (mub.*(Pressure-Pf)*sin(2*del+2*psi))./(sin(2*del+2*psi)+mub*(cos(2*del+2*psi)-cos(2*psi))); %shear stress due to internal friction (MPa) 
taup = 550+tfp; %shear stress (MPa)
%% overall strength due to anitgorite friction and semi-brittle strength (MPa)
tau_t = (ttl.^(-20)+tts.^(-20)+taup.^(-20)).^(-1/20);
%% Burdette and Hirth (2022) Low temperature Plasticity
Alp = exp(-0.624); %+/-.236 variables is flow law
t0 = 2.42e9; %+/- 0.09e9
qlp = 1.18; %+/- 0.09
Flp =86.3e3; % +/- 2.9e3
mulp = 35e9; %
sigma_Pa_lp = zeros(length(TK),1);
syms sigma;
for i = 1:length(TK)
sigma_Pa_lp(i) = vpasolve(e == Alp*(sigma/mulp)^2*exp(-1*Flp/(R*TK(i))*(1-(sigma/t0))^qlp), sigma, 100*10^6); %numerically solve for stress given strain rate (MPa)
end
tau_lp = sigma_Pa_lp/10^6/(3^0.5); %convert from axial symmetry to simple shear
%%
figure(1);
hold on
plot(tau_t(TC<650), -1*Depth(TC<650), tau_lp(TC<650), -1*Depth(TC<650), ttl(TC<400), -1*Depth(TC<400))
box on
xlim([0 500])
ylim([-80 -10])
legend('total antigorite friction', 'low temperature plasticity', 'lizardite friction' )
xlabel 'Shear Stress (MPa)'
ylabel 'Depth (km)'
%
%writematrix(ttl(TC<400),'Shumagin_EM.xlsx','Range','O1')
%data = [tau_weakest tau_lp]; %
%writematrix(data,'Shumagin_EM.xlsx','Range','P:Q')