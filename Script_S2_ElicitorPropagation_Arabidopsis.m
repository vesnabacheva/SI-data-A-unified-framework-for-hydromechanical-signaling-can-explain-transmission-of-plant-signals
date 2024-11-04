% Simulation of the advection-diffusion of Ricca factor in A.th. 
% This code was used to obtain the plot in Fig.4b
% September 2024, Vesna Bacheva

clear all
close all
clc


set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',1,...
          'defaultlinelinewidth',1,'defaultpatchlinewidth',1,'defaultaxesfontname','Times');
%% Paramaters 
D = 1e-10; % Elicitor diffusivity [m2 s-1] 
kx_ves = 2e-4; % Xylem vessel conductivity [mol m-1 Pa-1 s-1]
l_ves = 5e-6; % Xylem vessel diamater [m]
mu = 1.8*1e-5; % Molar volume of water [m3 mol-1]
L = 8*1e-3; % Leaf-to-leaf length [m]
T = 35; % Comsol simulation time [s]

%% Load Comsol data of gradient of xylem pressure 
NTP = load('File_F3_Xylem_Pressure_Gradient_Arabidopsis_E0_0.mat'); % Non-transpiring plant
grad_x_NT = NTP.grad_x_NT;

TP = load('File_F4_Xylem_Pressure_Gradient_Arabidopsis_E0_1.mat'); % Transpiring plant
grad_x_T = TP.grad_x_T;

time_length = size(grad_x_T,1); % Time points
x_length = size(grad_x_T,2); % Spatial points

dt  = T/time_length; % Time increments
dx  = L/x_length; % Time increments

t = 0:dt:T;
x = 0:dx:L;

%% Run Lagrangian simulation 
Np = 50; % number of elicitor particles;

%Line injection of Np elicitors at wound site, x = -L/2
%Transpiring plant
xpart_T(1:Np) = 0;
ypart_T = linspace(0,l_ves,Np);

%Location of elicitors
x_sim_T = zeros(length(t), Np);
y_sim_T = zeros(length(t), Np);

%Non-transpiring plant
xpart_NT(1:Np) = 0;
ypart_NT = linspace(0,l_ves,Np);

%Location of elicitors
x_sim_NT = zeros(length(t), Np);
y_sim_NT = zeros(length(t), Np);

for ti = 1:time_length
    % Xylem vessel velocity [m s-1]
    vx_T = -kx_ves.*grad_x_T(ti,:)*mu; % Transpiring plant
    vx_NT = -kx_ves.*grad_x_NT(ti,:)*mu; % Non-transpiring plant

    for i = 1 : Np
        % Transpiring plant
        index_T = (find(min(abs(x-xpart_T(i))) == abs(x-xpart_T(i))));
        Pe_T = vx_T(index_T)*l_ves/D;
        xpart_T(i) = vx_T(index_T).*dt + normrnd(xpart_T(i),sqrt(2.*dt.*D*(1+Pe_T^2/48)));
        ypart_T(i) = normrnd(ypart_T(i),sqrt(2.*dt.*D*(1+Pe_T^2/48)));

        % Non-transpiring plant
        index_NT = (find(min(abs(x-xpart_NT(i))) == abs(x-xpart_NT(i))));
        Pe_NT = vx_NT(index_NT)*l_ves/D;
        xpart_NT(i) = vx_NT(index_NT).*dt +normrnd(xpart_NT(i),sqrt(2.*dt.*D*(1+Pe_NT^2/48)));
        ypart_NT(i) = normrnd(ypart_NT(i),sqrt(2.*dt.*D*(1+Pe_NT^2/48)));
    end


x_sim_T(ti,:)= xpart_T;
y_sim_T(ti,:)= ypart_T;

x_sim_NT(ti,:)= xpart_NT;
y_sim_NT(ti,:)= ypart_NT;

end

%% Plot predicted displacement 
figure(1)
plot(t, mean(x_sim_T,2).*1e2, '-k')
hold on
plot(t, mean(x_sim_NT,2).*1e2, '--k')
xlabel('Time from wound [s]')
ylabel('Distance [cm]')
xlim([0 30])
ylim([0 0.75])
legend('Transpiring plant', 'Non-transpiring plant')
set(gcf,'color','w');
hold on