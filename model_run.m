%% Load and calculate WT ER geometry

% Load mesh of WT ER network from simulation output file
dirname = '/data/proj/ERCaSims/results/refill/diamond3degNF/start0/';
filename = 'diamondR10_nf_nuccon91.9.0.mesh.txt';
filename = [dirname filename];
MSH = MeshObj(filename);

a = 0.05; % ER tubule radius (μm)

% Compute total ER volume assuming all segments are cylindrical
totvol_WT = sum(MSH.len) * pi * a^2;

% Identify nuclear region (resvind == 1 marks nuclear sheet)
isnuc = (MSH.resvind == 1);

% Compute total ER surface area:
% - For tubules: side area of cylinders (excluding nuclear sheet)
% - For nuclear sheet: approximate as a sphere with radius 5 μm
totA_WT = sum(MSH.len(~isnuc)) * 2 * pi * a + 4 * pi * 5^2;

% Compute volume of peripheral (non-nuclear) ER
periphvol_WT = sum(MSH.len(~isnuc)) * pi * a^2;

%% Load and calculate RTN3OE ER geometry

% Load mesh of RTN3OE ER with vesiculated structure
dirname = '~/UCSD/data/ERCaSims/results/refill/HCPbub/deg6_bubR0pt6tubeRpt018lenpt3_10sheet_buffscl/';
filename = 'HCP_nf_bub_10sheet_buffscl.2.0.mesh.txt';
filename = [dirname filename];
MSH = MeshObj(filename);

a = 0.018; % ER tubule radius in RTN3OE (μm)

% Compute total ER volume (same formula)
totvol_RTN3OE = sum(MSH.len) * pi * a^2;

% Identify nuclear and vesicle regions
isnuc = (MSH.resvind == 1);
isbub = (MSH.resvind > 1);
nbub = nnz(isbub);  % Number of vesicles

% Estimate vesicle radius assuming spherical geometry
bubrad = (MSH.len(find(isbub,1)) / (4/3 * pi) * pi * a^2)^(1/3);

% Compute total surface area:
% - Tubules (excluding nuclear and vesicle)
% - Vesicles as spheres
% - One nuclear sheet approximated as 5 μm radius sphere
totA_RTN3OE = sum(MSH.len(~isbub & ~isnuc)) * 2 * pi * a + nbub * 4 * pi * bubrad^2 + 4 * pi * 5^2;

% Compute non-nuclear ER volume (tubules + vesicles)
bubvol_RTN3OE = sum(MSH.len(~isnuc)) * pi * a^2;

%% Compute volume scaling factor between WT and RTN3OE ER
volchange = totvol_RTN3OE / totvol_WT;

%% Setup Ca2+ model parameters

volchange = 1;  % Set to 1 to ignore geometric scaling in this run

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
params.k_soce = 6.5227;
%params.k_soce = 3.9583;

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
