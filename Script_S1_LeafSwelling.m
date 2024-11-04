% Simulation of leaf thickness increase following a wounding event
% This code was used to obtain Movie S4 and plots in Fig.3 b-c
% September 2024, Vesna Bacheva

clear all
close all
clc

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',1,...
          'defaultlinelinewidth',1,'defaultpatchlinewidth',1,'defaultaxesfontname','Times');

%%
model = createpde;

%paramaters
kappa_c = 1e-10; % Tissue poroelastic diffusivity [m2 s-1]
W = 200e-6; % Vein spacing  [m]
d_th = 200e-6; % Leaf  thickness [m]
tau_t = (W^2)/kappa_c; % Tissue poroelastic time scale [s]

tlist = [0:1:1200]; % Simulation time in [s]
%% Create geometry
rect1 = [3, 4, -W, W, W, -W, -d_th/2, -d_th/2, d_th/2, d_th/2]';
rect2 = [3, 4, -d_th/6, d_th/6, d_th/6, -d_th/6, -d_th/6, -d_th/6, d_th/6, d_th/6]';

%combine shapes
sf = 'rect1-rect2';
gm = [rect1, rect2];
ns = char('rect1','rect2');
ns = ns';

g = decsg(gm,sf,ns);

geometryFromEdges(model,g);
%pdegplot(model,'EdgeLabels','on','FaceLabels','on'); % Plot geometry

%% Generate the meash and specify model coefficients
m1 = generateMesh(model);
%pdeplot(m1) % Plot mesh 

% Specify the model 
specifyCoefficients(model,'m',0,'d',1,'c',kappa_c,'a',0,'f',0);

%% Set initial conditions
setInitialConditions(model,1); 

%Set Boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',[3 4 5 8],'u',0);
applyBoundaryCondition(model,'neumann','Edge',[1 2 6 7],'q',0,'g',0); 


%Solve the model on the time list tlist
result = solvepde(model,tlist);

%% plot the results
avg_Pt = zeros(1,size(tlist,2));% Leaf thickness increase in [%]

for ii=1:length(tlist)
  pdeplot(model,'XYData',result.NodalSolution(:,ii)*100);
  avg_Pt(ii) = mean(result.NodalSolution(:,ii))*100; 
  title(sprintf('Time = %.2f s', tlist(ii)));
  caxis([0 100]);
  colormap('jet')
  axis equal
  axis off
  set(gcf,'color','w');
  hcb=colorbar;
  title(hcb,'Relative dilation [%]')
  pause(0.1);
end


%% Plot resulting leaf thickness increase
figure(1)
plot(tlist./60, 100-avg_Pt, '-k');
xlim([-4 20]);
ylim([-10 110]);
xticks([-4 0 4 8 12 16 20]);
xlabel('t[min]');
ylabel('Change in normalized leaf thickness [%]')
title(sprintf('W = %.2f mm', W*1e3));
set(gcf,'color','w');
grid on
hold on

