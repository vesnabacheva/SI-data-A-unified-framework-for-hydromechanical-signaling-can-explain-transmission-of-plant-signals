% Simulation of the propagation of local tissue pressure following a wounding event
% This code was used to obtain the simulation shown in Movie S5
% September 2024, Vesna Bacheva

clear all
close all
clc

%% Parameters 
global kappa_t W P0
T = 200; % Time of interest [s]
W = 0.5*1e-3; % Vein spacing [m]
r_wound = 5e-6; % Radius of wounded cell [m]
kappa_t = 1e-10; % Tissue Poroelastic diffusivity [m2 s-1]
P0 = -0.2; % Initial tissue water potential [MPa]
pi_t = 0.75; % Tissue osmotic pressure [MPa]

%% Solve local tissue pressure
% Spatial and time discretization
nr = 50;     % Number of spatial points
nt = 50;     % Number of time points
r = linspace(r_wound, W, nr);  % Radial points from r_c to W
t = linspace(0, T, nt); % Time points from 0 to T
    
% Solve the PDE using pdepe
m = 1;  % Symmetry parameter for cylindrical coordinates
sol = pdepe(m, @radial_pde, @initial_condition, @boundary_conditions, r, t);

%% Plot the 2D radial solution (colormap)
figure(1);
set(gcf, 'Color', 'w');  % Set figure background to white
colormap('jet');
colormap(flipud(colormap));
hold on;
theta = linspace(0, 2*pi, 200); % Angular coordinates
[Theta, R] = meshgrid(theta, r); % Create polar grid
X = R .* cos(Theta); % Convert to Cartesian coordinates (x)
Y = R .* sin(Theta); % Convert to Cartesian coordinates (y)

% Plot a filled blue circle at the wound cell
x_circle = r_wound*1e6 * cos(theta); % X coordinates of the circle
y_circle = r_wound*1e6 * sin(theta); % Y coordinates of the circle
fill(x_circle, y_circle, [0, 0, 1], 'EdgeColor', 'none'); % Filled blue circle

for i = 1:length(t)
    % Get the solution at the selected time
    u_t = sol(i, :) + pi_t; % Tissue turgor pressure 
    
    % Plot the 2D radial solution
    Z = repmat(u_t', 1, length(theta)); % Replicate u(r) over all angles
    contourf(X*1e6, Y*1e6, Z, 80, 'LineStyle', 'none'); % Plot filled contour
    
    % Add title and formatting
    title(['Radial Diffusion at t = ' num2str(t(i)) ' s']);
    xlabel('x [\mum]');
    ylabel('y [\mum]');
    axis equal;
    colorbar;
    hcb = colorbar;
    title(hcb, 'Tissue turgor pressure [MPa]');
   
    axis off
    pause(0.5); % Pause to show the time evolution
end
hold off;

%% PDE functions
function [c, f, s] = radial_pde(r, t, u, dudr)
    global kappa_t
    c = 1;  
    f = kappa_t * dudr;  
    s = 0;  
end

% Initial condition function
function u0 = initial_condition(r)
global P0
    u0 = P0;
end

% Boundary conditions function
function [pl, ql, pr, qr] = boundary_conditions(~, ul, ~, ur, ~, r)
    % Boundary condition at r = r_c: u(r_c) = 0
    pl = ul;  
    ql = 0; 
    
    % Boundary condition at r = W: du/dr(W) = 0
    pr = 0;  
    qr = 1;   
end
