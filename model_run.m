% This script simulates Ca2+ dynamics with hardcoded parameters for two cases:
% - Wild Type (WT): set volchange = 1 and params.k_soce = 6.5227
% - RTN3OE: set volchange = 4.3808 and params.k_soce = 3.9583
% Ensure the matching k_soce value is uncommented when switching between cases

%% Setup Ca2+ model parameters

volchange = 1;  % Set to 1 to run simulation with the ER volume of WT

% Volume of ER and cytosol (μm^3)
params.VER = 212;
params.Vcyto = 4/3 * pi * (10^3 - 5^3) - params.VER;  % Cell minus nucleus minus ER

% IP3R Ca2+ release rate (μm^3/s)
params.v1 = 69.075;

% Passive leak rate (μm^3/s)
params.v2 = 0.00;

% SERCA pump rate (μM·μm^3/s)
params.v3 = 496;

% Cytoplasmic clearance rate (μM·μm^3/s)
params.v4 = 368;

% Half-max [Ca2+] for cytosolic clearance (μM)
params.k_pl = 0.106;

% SOCE entry rate (μm^3/s)
params.k_soce = 6.5227;  % for WT
%params.k_soce = 3.9583;  % for RTN3OE

% Extracellular Ca2+ concentration (μM)
params.c_cs = 200;

% Half-max constants for IP3R gating (μM)
params.k_er = 0.117;     % SERCA half-max [Ca2+] (μM)
params.d_ip3 = 0.1300;   % IP3 binding (μM)
params.d_act = 0.0823;   % Activation by Ca2+ (μM)
% inhibition by Ca, not using the simplified model, different IP3 binding
% strength depending on whether Ca is present in other sites
d2 = 1.049; d3 = 0.9434;
params.d_inh = d2*(params.IP3 + params.d_ip3)/(params.IP3 + d3); % μM
params.a = 0.5940; % h-gate kinetics (1/μM/s)

% IP3 concentration (μM)
params.IP3 = 0.4;

% Buffer parameters in ER lumen
params.k_D = 201.6;  % Dissociation constant (μM)
params.S = 2700 / volchange;  % Total buffer site conc. (μM)

params0 = params;  % Save a copy

%% Initial conditions
c_cyto0 = 0.06;      % Initial cytosolic free Ca2+ (μM)
h0 = 0.8;            % Initial h-inactivation variable
u_er0 = 170;         % Initial ER free Ca2+ (μM)
c_er0 = u_er0 + u_er0 * params.S / (u_er0 + params.k_D);  % Total ER Ca2+ (free + bound)

% Initial state vector: [cytosol Ca2+, h, ER free Ca2+, ER total Ca2+]
y0 = [c_cyto0; h0; u_er0; c_er0];

% Time span (s)
tspan = [0, 480];

% Solve ODE system
[t, y] = ode45(@(t, y) calcium_ode_cyto_V(t, y, params), tspan, y0);

%% Plot cytosolic and ER free Ca2+ over time
subplot(2,1,1);
plot(t/60, y(:,1), 'LineWidth', 2); 
plot_cleanup(gca, 'FontSize', 16)
ylabel('$c_\mathrm{cyto}\;\; \mu M$'); 
xlabel('Time (min)'); 
title('Cytosol Ca$^{2+}$');
legend off

subplot(2,1,2);
plot(t/60, y(:,3), 'LineWidth', 2); 
plot_cleanup(gca, 'FontSize', 16)
ylabel('$u_\mathrm{ER}\;\; \mu M$'); 
xlabel('Time (min)'); 
title('ER free Ca$^{2+}$');
legend off
