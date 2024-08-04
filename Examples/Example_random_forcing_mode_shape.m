% This script simulates an osillator chain with general damping and various
% modal decomposition methods are used to identify the modal parameters
% from the data. This script is used to generate comparitive results to the
% identified modes, poles, damping coeffs., and natural frequencies.
%
% Hewenxuan Li April 2023.
% Modifications:
% 2023-04-27: This code is dedicated to data generation purpose
% 2023-09-25: Updated code for noise level robustness
% 2024-03-10: updated the visualization for publication
% 2024-08-03: updated for release
%% Load Packages and Assign Paths
clear
close all
tic
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

%% MDOF system
RandomIC = 0;               % Using random initial condition?
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
% Additional dampers
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

%% Visualize the analytical modes
% Plot the complex normal modal
mode_indx = [1,2,3];
ComplexModePlot(modesn, mode_indx);

% Plot the complex normal modes in yyaxis form
mode_indx = [1,2,3];
ComplexModePlot_yy(modesn, mode_indx);

% Plot the Complex Normal modes in general complex form
PlotComplexModes(modesn)
% print('-dpng','-r600',strjoin({FigurePath,'Analytical_modes_PD.png'}, filesep))

%% Exogenous Loading function
close all
rr = 1;                           % Define the resampling rate (default is 1)
rseed = 0;                        % Specify random seed
rng(rseed);                       % Specify the seed for consistent results 
loading = 'Random';               % Specify the type of loads
timing.Tstart = 0;                % Specify the starting time of the simulation
timing.Tend = tend;               % Specify the termination time of the simulation
timing.dT = dt*rr;                % Samping period (after resampling)
finfo.A = 10;                     % Excitation amplitude 
finfo.omega = [];                 % Excitation frequency
finfo.windtype = 'Rectangular';   % Window type if windowing of the excitation is involved 
finfo.pwind = 0;                  % Period of windowing (default is 0, no windowing)
finfo.floc = 'beam';              % Excitation location (not in use in this example)
Tend = timing.Tend;               
RandomForcing = [];
for i = 1:n
[forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
RandomForcing = [RandomForcing, ppval(forcing, tspan)'];
end
disp('INSPECTION OF THE FORCING FUNCTION:')
DataInspection(ppval(forcing,tspan),tspan,fs)
set(gcf,'paperposition',[0,0,5,4])
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

% Add measurement noise
snr = inf;
seed = 0;
Y = awgn(Y,snr,'measured',seed);
dY = awgn(dY,snr,'measured',seed);

% DataInspection([real(q(:,[1,3])),real(dq(:,[1,3]))], tspan, fs)
% xlim([0 10])

DataInspection([real(Y(:,[1,3])),real(dY(:,[1,3]))], tspan, fs)
xlim([0 2])
set(gcf,'paperposition',[0,0,5,4])
%% SAVE SYSTEM INFORMATION AND 
SaveData = 0;
if SaveData == 1
% Specify the file name and directory path
filename = 'Simulation_Setup.txt';
SpecialCaseName = [];
DataName = ['Analytical_SNR_',num2str(round(snr))];
if ~isempty(SpecialCaseName)
    CaseName = strjoin({'MDOF','n',num2str(n),loading,DampType,'Ntrials',num2str(N),SpecialCaseName},'_');
else
    CaseName = strjoin({'MDOF','n',num2str(n),loading,DampType,'Ntrials',num2str(N)},'_');
end
directory = strjoin({DataPath,DataName,CaseName},filesep);
% Check if the directory exists
if ~exist(directory, 'dir')
    % If the directory does not exist, create it
    mkdir(directory);
end

% Save system setup
System.M = M; System.C = C; System.K = K;
System.RandomIC = RandomIC; System.dt = dt;
if RandomIC ~= 0
    System.Nic = N;
end
System.fs = fs; System.tspan = tspan; System.SystemOrder = n;
% Save data 
data.Y = Y; data.dY = dY;
save(strjoin({directory,'System.mat'},filesep),'System')
save(strjoin({directory,'data.mat'},filesep),'data')
disp(['System info Saved to:',strjoin({directory,'System.mat'},filesep),'System'])
disp(['Data saved to:',strjoin({directory,'data.mat'},filesep),'data'])
% Save the above setup figures
figure(1), print(strjoin({directory,strjoin({'Forcing',CaseName},'_')},filesep),'-r600','-dpng')
figure(2), print(strjoin({directory,strjoin({'Response',CaseName},'_')},filesep),'-r600','-dpng')
disp(['Plots saved to:',strjoin({directory})])

% NUMERICAL EXPERIMENT REPORT

% Create the full file path
full_path = fullfile(directory, filename);

% Open the file for writing
fid = fopen(full_path, 'w');

fprintf(fid,"====================== RESULTS FOR AMM DISCRETIZATION =================================\n");
fprintf(fid,"For Research Use: output-only modal analysis tool box v0.0\n");
fprintf(fid,"Copyright by Hewenxuan Li, hewenxuan.li@cornell.edu\n");
fprintf(fid,strcat('Report Generated on: ', string(datetime),'\r\n'));
fprintf(fid,"---------------------- MDOF SYSTEM SPECIFICATION ---------------------------------------------\n");
fprintf(fid,['System Order: ',num2str(n), '\n']);
fprintf(fid,['Damping Type: ', DampType, '\n']);

fprintf(fid,"-------------------------- STEPPING SETUP ---------------------------------------------\n");
fprintf(fid,['Sampling Rate (fs): ', num2str(fs), ' Hz\n']);
fprintf(fid,['Simulation Time (tend): ', num2str(tend), ' Secs\n']);
fprintf(fid,['Forcing Type: ', loading, '\n']);
% fprintf(fid,['Stepping Method: ', IntMethod, '\n']);
% fprintf(fid,['Weight of Forcing: ', ForceWeight, '\r\n']);

fprintf(fid,"---------------------- RESPONSE INFORMATION -------------------------------------\r\n");
fprintf(fid,'\n');
fprintf(fid,['The Size of the Configuration Data:',num2str(size(Y,1)),'-by-',num2str(num2str(size(Y,2))),'\n']);
% fprintf(fid,['The Manitude of Modal Participation (first 5 modes):',num2str(IW(2,1:5)),'\n']);
% fprintf(fid,['Modes Included: ', num2str(ModeIndices), '\n']);
fprintf(fid,"=======================================================================================");

fclose(fid);
disp('Report generated!')
disp(['Report Directory:',full_path])
end
toc
%% 
close all
mode_indx = [1,4,6,9];
Metric = 'rell2';
%% GCVD (Left EVP formulation)
% Conduct GCVD (left EVP formulation in the paper)
tic
[U, V, S, Lambda, W_gcvd] = gcvd(Y,dY,'forced');
T_gcvd = toc;
% Define the estimated modes
Phi_hat_gcvd = inv(W_gcvd');
% Sort the modes according to the frequencies
[Phi_gcvds, Lambda_gcvds] = modeSortMAC(Phi_hat_gcvd, S{2}/S{1}, modes(2:end-1,:), 'Uniquedisp');
% Plot preparation -> decomposing complex modes into real and imag parts
Phi_gcvd = ReImModes(Phi_gcvds, 'mdof_ff');
% Visualize the complex modes in 3D space
ComplexModeComp(Phi_gcvd, modes, mode_indx, 'default', ',\,gcvd')
% Visualize the MAC values
[mac_gcvd,~,~] = MAC(Phi_gcvds(1:10, :), modes(2:end-1, :));
set(gcf,'renderer','painters')
% Calculate the identification error
ephi_gcvd = ComplexModeError(modes, [zeros(1,n); Phi_gcvds; zeros(1,n)], Metric);

% Compare the poles
mps_gcvd = poles2mp(diag(Lambda_gcvds));
[~, Lambda_gcvds_full] = modeSortMAC(Phi_hat_gcvd, S{2}/S{1}, modes(2:end-1,:), 'Full');
PlotPoles({Lambdau,Lambda_gcvds_full},{'o','x'},{},[1,1],{'Truth','GCVD-A'})

% Adjust the estimated system poles
P1 = U*S{1};
% Visualize the adjusted pole estimates
Lambda_gcvd_bar = pinv(P1)*((dY-(F*B'))*W_gcvd); 
[~, Lambda_gcvd_bar_full] = modeSortMAC(inv(W_gcvd'), Lambda_gcvd_bar, modes(2:end-1,:), 'Full');
PlotPoles({Lambdau,Lambda_gcvds_full,Lambda_gcvd_bar_full},{'o','x','*'},{},[1,0.9,0.8],{'$\mathrm{Truth}$','$\Lambda_{gcvd}$','$\bar\Lambda_{gcvd}$'})
xlim([-0.7 0])

% Caculate the rectified poles
mps_gcvd_r = poles2mp(Lambda_gcvd_bar_full);
%% DMD
% Conduct DMD
tic
[mode_dmd, Lambda_dmd, b_dmd] = DMD(Y(1:end-1,:)', Y(2:end,:)',rank(Y));
T_dmd = toc;
% Ensure the esimated modes are complex-valued
mode_dmd = complexModeCheck(mode_dmd);
% Get the DMD continuous-time eigenvalue matrix
Lambda_dmd = diag(log(diag(Lambda_dmd))/dt);
% Sort the modes and obtain only the unique displacement modes
[Phi_dmds, Lambda_dmds] = modeSortMAC(mode_dmd, Lambda_dmd, modes(2:end-1,:),'Uniquedisp');
[mac_dmd,~,~] = MAC(Phi_dmds, modes(2:end-1, :));

Phi_dmd = ReImModes(Phi_dmds, 'mdof_ff');
ComplexModeComp(Phi_dmd, modes, mode_indx, 'default',',\,dmd')
% Calculate the identification error
ephi_dmd = ComplexModeError(modes, [zeros(1,n); Phi_dmds; zeros(1,n)], Metric);

% Compare the poles
mps_dmd = poles2mp(diag(Lambda_dmds));
[~, Lambda_dmds_full] = modeSortMAC(mode_dmd, Lambda_dmd, modes(2:end-1,:),'Full');
PlotPoles({Lambdau,Lambda_dmds_full},{'o','x'},{},[1,1],{'Truth','DMD'})
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
tic
[ERA_Result] = ERA(Data,Sampling_Frequency,Hankel_nColumns,Hankel_nRows,1,Cutoff,Shift,EMAC_option);
T_era = toc;
ERA_Mode_Shapes = ERA_Result.Parameters.ModeShape;
ERA_Poles = diag(ERA_Result.Parameters.Poles);
% Ensure the esimated modes are complex-valued
ERA_Mode_Shapes = complexModeCheck(ERA_Mode_Shapes);
% Sort the modes and poles using MAC method
[Phi_eras, Lambda_eras] = ModeSortFreq(ERA_Mode_Shapes, ERA_Poles, 'Uniquedisp');
[mac_era,~,~] = MAC(Phi_eras, modes(2:end-1, :));

mps_era = poles2mp(diag(Lambda_eras));
% Plot preparation -> decomposing complex modes into real and imag parts
Phi_era = ReImModes(Phi_eras, 'mdof_ff');
% Calculate the identification error
ephi_era = ComplexModeError(modes, [zeros(1,n); Phi_eras(:,1:n); zeros(1,n)], Metric);
% Compare the poles
[~, Lambda_era_full] = modeSortMAC(ERA_Mode_Shapes, ERA_Poles, modes(2:end-1,:), 'Full');
PlotPoles({Lambdau,Lambda_era_full},{'o','x'},{},[1,1],{'Truth','ERA'})
% Visualize the complex modes in 3D space
ComplexModeComp(Phi_era, modes, mode_indx, 'default', ',\,era')
set(gcf,'paperposition',[0,0,5,3])
%% ITD
% Conduct ITD
tic
[Lambda_itd, mode_itd] = ITD({Y});
T_itd = toc;
% Get the ITD continuous-time eigenvalue matrix
Lambda_itd = diag(log(diag(Lambda_itd))/dt);
% Sort the modes and obtain only the unique displacement modes
[Phi_itds, Lambda_itds] = modeSortMAC(mode_itd(1:20,:), Lambda_itd, modes(2:end-1,:),'Uniquedisp');

% Isolate the displacement modes
% Phi_itds = Phi_itds(n+1:end,:);
Phi_a = Phiu(1:n,:);
[mac_itd,~,~] = MAC(Phi_itds, modes(2:end-1,:));
Phi_itd = ReImModes(Phi_itds, 'mdof_ff');
ComplexModeComp(Phi_itd, modes, mode_indx, 'default',',\,itd')
% Calculate the identification error
ephi_itd = ComplexModeError(modes, [zeros(1,n); Phi_itds(:,1:n); zeros(1,n)], Metric);

% Compare the poles (this shows that the poles identified by ITD are not necessarily the physical poles of the system)
mps_itd = poles2mp(diag(Lambda_itds));
mps_itd = mps_itd(:,1:n);
[~, Lambda_itds_full] = modeSortMAC(mode_itd(1:20,:), Lambda_itd, modes(2:end-1,:), 'Full');
PlotPoles({Lambdau,Lambda_itds_full},{'o','x'},{},[1,1],{'Truth','ITD'})
% Zoom to see the identified poles that are close to the true poles
PlotPoles({Lambdau,Lambda_itds_full},{'o','x'},{},[1,1],{'Truth','ITD'})
xlim([-0.7 0.4])
ylim([-4 4])
%% Visualize the Mode Shape Comparison
% Check file system and create directory for plot export
ExamplePath = strjoin({FigurePath,'ModalParameters'},filesep);
datetime.setDefaultFormats('default','MMM-dd-yyyy');
[ret, HostName] = system('hostname');
% Filter out unwanted hostname elements (e.g., return)
HostNameCheck = split(HostName);
if length(HostNameCheck) > 1
    HostName = HostNameCheck{1};
end

TrialName = strjoin({char(datetime),HostName},'-');
TrialPath = strjoin({ExamplePath, TrialName}, filesep);
NewTrial = 'False';

% Examine if there exists the trial path
% Check if the directory exists
counter = 1;
if exist(TrialPath, 'dir') ~= 7 % If there does not exist a trial path
    % If it doesn't exist, create the directory
    mkdir(TrialPath);
    disp(['Directory created: ', TrialPath]);
else % If it exists
    if strcmpi(NewTrial, 'true') % If a new trial is wanted, create a new one by lumping new trial numbers after the path
        disp(['Directory already exists: ', TrialPath]);
        tempNewPath = strjoin([TrialPath, string(counter)],'-');
        while exist(tempNewPath, 'dir') == 7
            disp(['Directory already exists: ', tempNewPath]);
            counter = counter + 1;
            tempNewPath = strjoin([TrialPath, string(counter)],'-');
        end
        mkdir(tempNewPath);
        disp(['Directory already exists: ', tempNewPath]);
    elseif strcmpi(NewTrial, 'false') % If a new trial is not desired, use the existing one instead
        disp(['Directory already exists: ', TrialPath]);
        disp('Will use the exising one for figure storage')
    end
end

markers = {'^','v','d','s'};
mode_indx = [1,2,4];
colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
alphas = [1,1,1,1,1];
legends = {'1','2','3','4','5'};
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'default', {'GCVD','DMD','ERA','ITD'}, markers, colors)


mode_indx = [1];
colors = {'#0072BD','#0072BD','#0072BD','#0072BD','#0072BD'};
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'default', {'GCVD','DMD','ERA','ITD'}, markers, colors)
set(gcf,'papersize',[5 2.5])
set(gcf,'paperposition', [0,0,5,2.5])
legend('Location','NorthEast')
PlotName = 'Fig_Mode1_orthographic_comp';
mysavefig(TrialPath,PlotName)

mode_indx = [1];
colors = {'#0072BD','#0072BD','#0072BD','#0072BD','#0072BD'};
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'complexplane', {'GCVD','DMD','ERA','ITD'}, markers, colors)
legend('off')
set(gcf,'papersize',[3 3])
set(gcf,'paperposition', [0,0,3,3])
PlotName = 'Fig_Mode1_complexplane_comp';
mysavefig(TrialPath,PlotName)

mode_indx = [2];
colors = {'#D95319','#D95319','#D95319','#D95319','#D95319'};
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'default', {'GCVD','DMD','ERA','ITD'}, markers, colors)
set(gcf,'papersize',[5 2.5])
set(gcf,'paperposition', [0,0,5,2.5])
legend('Location','NorthEast')
PlotName = 'Fig_Mode2_orthographic_comp';
mysavefig(TrialPath,PlotName)
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'complexplane', {'GCVD','DMD','ERA','ITD'}, markers, colors)
legend('off')
set(gcf,'papersize',[3 3])
set(gcf,'paperposition', [0,0,3,3])
PlotName = 'Fig_Mode2_complexplane_comp';
mysavefig(TrialPath,PlotName)

mode_indx = [4];
colors = {'#EDB120','#EDB120','#EDB120','#EDB120','#EDB120'};
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'default', {'GCVD','DMD','ERA','ITD'}, markers, colors)
set(gcf,'papersize',[5 2.5])
set(gcf,'paperposition', [0,0,5,2.5])
legend('Location','NorthEast')
PlotName = 'Fig_Mode4_orthographic_comp';
mysavefig(TrialPath,PlotName)
ComplexModeComp_({Phi_gcvd,Phi_dmd,Phi_era,Phi_itd}, modes, mode_indx, 'complexplane', {'GCVD','DMD','ERA','ITD'}, markers, colors)
legend('off')
set(gcf,'papersize',[3 3])
set(gcf,'paperposition', [0,0,3,3])
PlotName = 'Fig_Mode4_complexplane_comp';
mysavefig(TrialPath,PlotName)

%% Pole Visualization
colors = {'#000000','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
PlotPoles({Lambdau,Lambda_gcvds_full,Lambda_dmds_full,Lambda_era_full,Lambda_itds_full},...
{'o','^','v','d','s'},...
colors,...
[0.7,1,0.8,0.8,0.8],...
{'Truth','GCVD','DMD','ERA','ITD'})
xlim([-1 0.3])
ylim([-4, 4])
set(gcf,'papersize',[3,3])
set(gcf,'paperposition',[0,0,3,3])
PlotName = 'Fig_pole_complexplane_comp_snr_inf';
mysavefig(TrialPath,PlotName)
%% Calculate  and visualize the mode shape errors
figure(2), clf
Metric = 'rell2';
PlotModeError({ephi_gcvd,ephi_dmd,ephi_era,ephi_itd},Metric, ...
    {'o','o','o','o'}, ...             % Define the markers
    {}, ...                                % Use default colors
    [0.5, 0.5, 0.5,0.5], ...          % Define alphas
    {'GCVD','DMD','ERA','ITD'})      % Define Legends
set(gca,'yscale','log')
%% Modal asurance comparison
figure(3),clf
plotMAC({mac_gcvd, mac_dmd, mac_era, mac_itd}, ...
    {'o-','o-','o-','o-'}, ...             % Define the markers
    {}, ...                                % Use default colors
    [0.5, 0.5, 0.5,0.5], ...          % Define alphas
    {'GCVD','DMD','ERA','ITD'})      % Define Legends
%% System modal damping ratio and natural frequency estimation error
plotMps({mps_gcvd, mps_dmd, mps_era, mps_itd}, mps, ...
    {'o-','o-','o-','o-'}, ...             % Define the markers
    {}, ...                                % Use default colors
    [0.5, 0.5, 0.5,0.5], ...          % Define alphas
    {'GCVD','DMD','ERA','ITD'})      % Define Legends
subplot(211), set(gca,'yscale','log')
%% With rectification for cvd and gcvd (equivalent to a EMA not OMA)
plotMps({mps_gcvd_r(:,1:2:end), mps_dmd, mps_era, mps_itd}, mps, ...
    {'o-','o-','o-','o-'}, ...             % Define the markers
    {}, ...                                % Use default colors
    [0.5, 0.5, 0.5,0.5], ...          % Define alphas
    {'GCVD','DMD','ERA','ITD'})      % Define Legends
subplot(211), set(gca,'yscale','log')
