% This example illustrates various data-driven modal analysis methods on
% estimating modal parameters of generally-damped systems 
% 
% Created by Hewenxuan Li April 2023.
% Contributors:
% Hewenxuan Li - dynamic simulation, gcvd, itd, dmd, and visualization
% Dalton Stein - era implementation
%
% Reference: Characteristic Value Decompositoin: A Unifying Paradigm for
% Data-Driven Modal Analysis, 2024
%
% Modifications:
% 2023-04-27: This code is created to obtain data for cvd development
% 2023-9-17: incorporated multiple methods, like itd, dmd, and era 
% 2024-08-03: updated for release
%% Load Packages and Assign Paths
clear
clc
close all
ProjName = detectParentFolder();
folders = split(pwd,filesep);
for i = 1:length(folders)
    if isequal(folders{i}, ProjName)
        flag = i;
        break
    end
end
ProjectPath = strjoin(folders(1:flag),filesep);
FigurePath = strjoin({ProjectPath,'Figures','MDOF'}, filesep);
DataPath = strjoin({ProjectPath,'Data','MDOF'}, filesep);
CodePath = strjoin({ProjectPath,'Code'}, filesep);
addpath(CodePath)
addpath(strjoin({ProjectPath, 'OS'},filesep))
addpath(strjoin({CodePath, 'Visualization'},filesep))
addpath(strjoin({CodePath, 'Metrics'},filesep))
% Plot Styles
checkdir(FigurePath);
[Markers, LineStyles, MyColors] = mystyle();
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultTextInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% MDOF system (10 DOF system, identical to the one used in the cvd paper)
fprintf('--------------- Define MDOF System -----------------\n')
RandomIC = 0;               % Using random initial condition or not
dt = 0.005;                 % Sampling time
fs = 1/dt;                  % Sampling Freq
tend = 150;                 % Simulation Time
tspan = 0:dt:tend;          % Discrete time vector
n = 10;                     % DOFs

% System Parameters
alpha = 1;                  % Parameters used for proportional damping
beta = 0;              
m = ones(n,1);              % Point mass
c = 0.3*ones(n+1,1);        % Damping coefficient
k = 4*ones(n+1,1);          % Stiffness
% Additional dampers that grants general damping
Ic = zeros(n+1,1);
index_Ic = [1 2];
Ic([1 2]) = 1;
ci = 0.3*ones(n+1,1);

% Create the system matrices M, C, and K from the fixed-fixed oscillator
% chain model
[M, C, K] = mdof_ff(n, m, c, k, index_Ic, ci);

% Check system damping type to be proportional or not
DampType = CheckDamping(M, C, K, n);

% Configuration Space Formulation - Generalized Eigenvalue Problem
[Phi, Lambda] = eig(inv(M)*K);
mps_ud = poles2mp(Lambda, 'Unique');
fprintf('\n')
disp(['Natural Frequencies: ',num2str(mps_ud(2,:))])

if isequal(lower(DampType),'md')
    zeta = [0.01, 0.1];
    Cm = 2*diag(zeta)*sqrt(Lambda);
    C = Phi*Cm*Phi;
elseif isequal(lower(DampType),'pd')
    alpha = 0.05;
    beta = 0.05;
    C = alpha*M + beta*K;
end

% Define the mode shape matrix (not normalized)
Phi_ud = [1 1; 1 -1];

%% Mirovitch/Robin Formulation - Eigenvalue problem
fprintf('------------ Obtain the ground truth modes --------------\n')
% Define system matrix
A = [zeros(n,n), eye(n,n); -inv(M)*K, -inv(M)*C];
B = [zeros(n,n); inv(M)];
CC = eye(size(A));    % Full observability
D = zeros(size(B));   % No control variable

[Phiu, Lambdau] = eig(A);
[Phiu, Lambdau] = ModeSortFreq(Phiu, Lambdau, 'Full');

mps = poles2mp(Lambdau, 'Unique');

[modes, poles] = ModeSortFreq(Phiu, Lambdau, 'Uniquedisp');
disp(['Damped Natural Frequencies: ',num2str(mps(2,:).*sqrt(1-mps(3,:).^2))])
disp(['Modal Damping Ratios: ',num2str(mps(3,:))])
modes = [zeros(1,size(modes,2));
         modes;
         zeros(1,size(modes,2))];
modesn = zeros(size(modes));
for i = 1:size(modes,2)
    modesn(:,i) = modes(:,i)/norm(modes(:,i));
end

%% Plot the complex normal modal
fprintf('----------- Visualize the Complex Mode Shapes --------------\n')
%% Plot the Complex Normal modes in general complex form
PlotComplexModes(modesn)
%% Exogenous Loading function
close all
rr = 1;
rng(1);                           % Specify the seed for consistent results 
loading = 'Random';
timing.Tstart = 0;
timing.Tend = tend;
timing.dT = dt*rr;
finfo.A = 10;
finfo.omega = [];
finfo.windtype = 'Rectangular';
finfo.pwind = 0;
finfo.floc = 'beam';
Tend = timing.Tend;
RandomForcing = [];
for i = 1:n
[forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
RandomForcing = [RandomForcing, ppval(forcing, tspan)'];
end
disp('==================== Analysis Parameters Updated ====================')
%% Forced oscillation
sys = ss(A, B, CC, D);
weight = zeros(1, n);     % Initialize forcing weight as zeros (no forcing)
If = [1 2 3 4 5 6 7 8 9 10];                 % Forcing index that stipulate on which mass the force is applied
weight(If) = 1;           % Assign the weights so that If-indexed mass applies the force
f = ppval(forcing,tspan)';% Forcing vector
% Define the forcing matrix F
if isequal(loading, 'Random')
    F = RandomForcing;        % If a random forcing, use the random inputs directly 
else
    F = repmat(f, 1, n);      % If no a random forcing, Matrix form of the forcing vector is diretly replicated
end
weight = repmat(weight, size(F,1), 1); % Forcing weight
F = F.*weight;
mode_weights = ones(n,1);
% mode_weights = randn(n, 1);
if RandomIC == 0
    N = 1;
    y0 = zeros(2*n,1);
    for i = 1:n
        y0 = y0 + mode_weights(i)*(Phiu(:,2*i-1) + Phiu(:,2*i)); % Define the initial displacement as a linear combination of the
        % modes (vectorize it later!)
    end
    y0(n+1:2*n,:) = 0; % Force the initial velocity to be zeros

    Y = lsim(sys,F,tspan,y0); % Simulate the system response
    
    % obtain the analytical velocity vector
    dY = A*Y' + B*F';
    dY = dY';
    
elseif RandomIC == 1 % (This part needs development)
    N = 5;
    y_temp = cell(N,1);
    dy_temp = cell(N,1);
    for i = 1:N
        alpha = randn(1,1); beta = randn(1,1);
        y0 = alpha*Phi_ud(:,1) + beta*Phi_ud(:,2);
        y0 = [y0; 0; 0]; 
        y_temp{i} = lsim(sys,F,tspan,y0);
        % obtain the analytical velocity vector
        dy = A*y_temp{i}' + B*F';
        dy = dy';
        dy_temp{i} = dy;
    end
    
    Y = cell2mat(y_temp);
    dY = cell2mat(dy_temp);
end
% Obtain the modal motion
q = inv(Phiu)*Y';
dq = inv(Phiu)*dY';
q = q';
dq = dq';

% DataInspection([real(q(:,[1,3])),real(dq(:,[1,3]))], tspan, fs)
% xlim([0 10])

DataInspection([real(Y(:,[1,3])),real(dY(:,[1,3]))], tspan, fs)
xlim([0 2])
set(gcf,'paperposition',[0,0,5,4])
%% Compare the mode decomposition methods
fprintf('----------- Compare the various methods ------------\n')
close all
mode_indx = [1,2,3,4];
Metric = 'rell2';
%% CVD (via GCVD and Left EVP formulation)
% Conduct GCVD (left EVP formulation in the paper)
[U, V, S, Lambda, W_gcvd] = gcvd(Y,dY,'forced');
% Sort the modes according to the frequencies
[Phi_gcvds, Lambda_gcvds] = ModeSortFreq(inv(W_gcvd'), S{2}/S{1}, 'Uniquedisp');

% Plot preparation -> decomposing complex modes into real and imag parts
Phi_gcvd = ReImModes(Phi_gcvds, 'mdof_ff');

% Calculate the identification error
ephi_gcvd = ComplexModeError(modes, [zeros(1,n); Phi_gcvds; zeros(1,n)], Metric);

% Compare the poles
mps_gcvd = poles2mp(diag(Lambda_gcvds));
[~, Lambda_gcvds_full] = ModeSortFreq(inv(W_gcvd'), S{2}/S{1}, 'Full');

% Adjust the estimated system poles
P1 = U*S{1};
% Visualize the adjusted pole estimates
Lambda_gcvd_bar = pinv(P1)*((dY-(F*B'))*W_gcvd); 
[~, Lambda_gcvd_bar_full] = ModeSortFreq(inv(W_gcvd'), Lambda_gcvd_bar, 'Full');
fprintf('CVD done ...\n')
%% DMD
% Conduct DMD
[mode_dmd, Lambda_dmd, b_dmd] = DMD(Y(1:end-1,:)', Y(2:end,:)',rank(Y));
% Get the DMD continuous-time eigenvalue matrix
Lambda_dmd = diag(log(diag(Lambda_dmd))/dt);
% Sort the modes and obtain only the unique displacement modes
[Phi_dmds, Lambda_dmds] = ModeSortFreq(mode_dmd, Lambda_dmd, 'Uniquedisp');

Phi_dmd = ReImModes(Phi_dmds, 'mdof_ff');
% Calculate the identification error
ephi_dmd = ComplexModeError(modes, [zeros(1,n); Phi_dmds; zeros(1,n)], Metric);

% Compare the poles
mps_dmd = poles2mp(diag(Lambda_dmds));
[~, Lambda_dmds_full] = ModeSortFreq(mode_dmd, Lambda_dmd, 'Full');
fprintf('DMD done ...\n')

%% ERA
% ERA Setup
Sampling_Frequency = 1/dt;
Maximum_Lag = 575; % Lags for Cross-Correlation using in NexTT (Not Used Currently)
Hankel_nColumns = 20000; % Columns in Hankel Matrix
Hankel_nRows =  20;   % Rows in Hankel
Cutoff = 20; % Cuttoff
Shift = 10; 
Data = Y(:,1:20);
Reference_Channels = 5;
EMAC_option = 0;
% Perform ERA
[ERA_Result] = ERA(Data,Sampling_Frequency,Hankel_nColumns,Hankel_nRows,1,Cutoff,Shift,EMAC_option);
ERA_Mode_Shapes = ERA_Result.Parameters.ModeShape;
ERA_Poles = diag(ERA_Result.Parameters.Poles);
[Phi_eras, Lambda_eras] = ModeSortFreq(ERA_Mode_Shapes, ERA_Poles, 'Uniquedisp');
% Plot preparation -> decomposing complex modes into real and imag parts
Phi_era = ReImModes(Phi_eras, 'mdof_ff');
ephi_era = ComplexModeError(modes, [zeros(1,n); Phi_eras(:,1:n); zeros(1,n)], Metric);
% Compare the poles
Lambda_era_full = sort(diag(ERA_Poles),'ascend');
fprintf('ERA done ...\n')

%% ITD
[Lambda_itd, mode_itd] = ITD({Y});
% Get the ITD continuous-time eigenvalue matrix
Lambda_itd = diag(log(diag(Lambda_itd))/dt);
% Sort the modes and obtain only the unique displacement modes
[Phi_itds, Lambda_itds] = ModeSortFreq(mode_itd, Lambda_itd, 'Uniquedisp');
% Isolate the displacement modes
Phi_itds = Phi_itds(n+1:end,:);
Phi_itd = ReImModes(Phi_itds, 'mdof_ff');
% Calculate the identification error
ephi_itd = ComplexModeError(modes, [zeros(1,n); Phi_itds(:,1:n); zeros(1,n)], Metric);
% Compare the poles (this shows that the poles identified by ITD are not necessarily the physical poles of the system)
mps_itd = poles2mp(diag(Lambda_itds));
[~, Lambda_itds_full] = ModeSortFreq(mode_itd, Lambda_itd, 'Full');
fprintf('ITD done ...\n')

%% Calculate the errors
PlotModeError({ephi_gcvd,ephi_dmd,ephi_itd,ephi_era},Metric,{'o','o','o','o'},{},[0.5,0.5,0.5, 0.5],{'GCVD','DMD','ITD','ERA'})
PlotPoles({Lambdau,Lambda_gcvd_bar_full, Lambda_dmds_full, Lambda_era_full, Lambda_itds_full},{'o','x','*','^','s'},{},[1,1,1,1,1],{'Truth','GCVD','DMD','ITD','ERA'})
xlim([-0.7 0.5])
ylim([-4 4])