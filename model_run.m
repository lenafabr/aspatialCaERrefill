% This script simulates Ca2+ dynamics with hardcoded parameters for two cases:
% - Wild Type (WT): set volchange = 1 and params.k_soce = 6.5227
% - RTN3OE: set volchange = 4.3808 and params.k_soce = 3.9583
% Ensure the matching k_soce value is uncommented when switching between cases

%% Setup Ca2+ model parameters

volchange = 1;  % Set to 1 to run simulation with the ER volume of WT
% volchange = 4.3808;  % Set to 4.3808 which is the ratio of ER volume of RTN3OE and WT 

% Volume of ER and cytosol (μm^3)
params.VER = 211 * volchange;
params.Vcyto = 4/3 * pi * (10^3 - 5^3) - params.VER;  % Cell minus nucleus minus ER

% IP3R Ca2+ release rate (μm^3/s)
params.v1 = 33;

% Passive leak rate (μm^3/s)
params.v2 = 0.00;

% SERCA pump rate (μM·μm^3/s)
params.v3 = 310;

% Cytoplasmic clearance rate (μM·μm^3/s)
params.v4 = 230;

% Half-max [Ca2+] for cytosolic clearance (μM)
params.k_pl = 0.114;

% SOCE entry rate (μm^3/s)
params.k_soce = 6.5227;  % for WT
%params.k_soce = 3.9583;  % for RTN3OE

% Extracellular Ca2+ concentration (μM)
params.c_cs = 200;

% Half-max constants for IP3R gating (μM)
params.k_er = 0.117;     % SERCA half-max [Ca2+] (μM)
params.d_ip3 = 0.1300;   % IP3 binding (μM)
params.d_act = 0.0823;   % Activation by Ca2+ (μM)
params.d_inh = 0.4139;   % Inhibition by Ca2+ (μM)
params.a = 0.2000 * 1.65; % h-gate kinetics (1/μM/s)

% IP3 concentration (μM)
params.IP3 = 0.5;

% Buffer parameters in ER lumen
params.k_D = 201.6;  % Dissociation constant (μM)
params.S = 2710 / volchange;  % Total buffer site conc. (μM)

params0 = params;  % Save a copy

%% Initial conditions
c_cyto0 = 0.06;      % Initial cytosolic free Ca2+ (μM)
h0 = 0.8;            % Initial h-inactivation variable
u_er0 = 182;         % Initial ER free Ca2+ (μM)
c_er0 = u_er0 + u_er0 * params.S / (u_er0 + params.k_D);  % Total ER Ca2+ (free + bound)

% Initial state vector: [cytosol Ca2+, h, ER free Ca2+, ER total Ca2+]
y0 = [c_cyto0; h0; u_er0; c_er0];

% Time span (s)
tspan = [0, 1000];

% Solve ODE system
[t, y] = ode45(@(t, y) calcium_ode_cyto_V(t, y, params), tspan, y0);

%% Plot cytosolic and ER free Ca2+ over time
subplot(2,1,1);
plot(t, y(:,1), 'LineWidth', 2); 
ylabel('c_{cyto} μM'); 
title('Cytosol Free Ca^{2+}');
plot_cleanup(gca, 'FontSize', 16)
legend off

subplot(2,1,2);
plot(t, y(:,3), 'LineWidth', 2); 
ylabel('u_{ER} μM'); 
xlabel('Time'); 
title('ER Free Ca^{2+}');
plot_cleanup(gca, 'FontSize', 16)
legend off
