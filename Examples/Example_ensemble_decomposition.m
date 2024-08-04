% This script simulates an osillator chain with general damping. 
% It conduct multiple trials of numerical experiments with distinct:
% 1. Forcing amplitudes
% 2. Signal-to-noise ratios
% 3. Sampling frequencies and
% 4. Experiment time periods
% This code conduct the experiments and compare the results over different
% numbers of trials N = 2^{0, 1, 2, 3, 4}. The per-mode error plots are
% shown in an overlay of line and scatter plots and a summary of a
% mode-wise averaged plot is shown in the very end of the code.
% Hewenxuan Li, 2023.
% Modifications:
% v0.1 - 2023-04-27: This code is dedicated to data generation purpose
% v1.0 - 2023-10-05: Updated for noise robustness exploration
% v2.0 - 2023-11-30: Updated the comparison between GCVD and ITD using
% Ensemble quotient matrices
% v2.1 - 2024-05-23: Added notes and cleaned up the research code
%% Load Packages and Assign Paths
clc
clear
close all
disp('=========== CVD: Ensemble Decomposition Example ==============')
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
disp('------------ MDOF Definition ---------------')
RandomIC = 0;               % Using random initial condition?
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
disp('Compute undamped system poles from EVP:')
disp(['Natural Frequencies: ',num2str(mps_ud(2,:))])
fprintf('\n')
disp('System damping type:')
if isequal(lower(DampType),'md')
    zeta = [0.01, 0.1];
    Cm = 2*diag(zeta)*sqrt(Lambda);
    C = Phi*Cm*Phi;
    disp('The system is modally damped.')
elseif isequal(lower(DampType),'pd')
    alpha = 0.05;
    beta = 0.05;
    C = alpha*M + beta*K;
    disp('The system is proportionally damped.')
elseif isequal(lower(DampType),'np')
    disp('The system is generally damped.')
end

% Define the mode shape matrix (not normalized)
Phi_ud = [1 1; 1 -1];

%% Mirovitch/Robin Formulation - Eigenvalue problem
fprintf('\n')
disp('Constructing the state space model.')
% Define system matrix
A = [zeros(n,n), eye(n,n); -inv(M)*K, -inv(M)*C];
B = [zeros(n,n); inv(M)];
CC = eye(size(A));    % Full observability
D = zeros(size(B));   % No control variable

[Phiu, Lambdau] = eig(A);
[Phiu, Lambdau] = ModeSortFreq(Phiu, Lambdau, 'Full');

mps = poles2mp(Lambdau, 'Unique');

fprintf('\n')
disp('System damping information:')
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

%% Obtain Candidate Experiments
fprintf('\n')
disp('--------------- Ensemble Decomposition Started ------------------')
% Initialize the data cells
max_N_exp = 5;
n_exp = [1,2,4,8,16];
% ------------------------------------------
% Print the total number of experiments
fprintf('Total number of experiments -> %d;\n', max_N_exp);

% Loop through each group and print the number of experiments
for i = 1:max_N_exp
    fprintf('Group %d has %d experiment%s', i, n_exp(i), (n_exp(i) > 1) * 's' + ''); % Add 's' if more than one experiment
    fprintf('.\n');
end
% ------------------------------------------

% Mode error cell
ephi_cell_cvd = cell(max_N_exp,1);
ephi_cell_itd = cell(max_N_exp,1);
% MAC cell
mac_cell_cvd = cell(max_N_exp,1);
mac_cell_itd = cell(max_N_exp,1);
% MPS cell
mps_cell_cvd = cell(max_N_exp,1);
mps_cell_itd = cell(max_N_exp,1);

for k = 1:max_N_exp
    N_exp = n_exp(k);    % Number of experiments
    fprintf('Conducting group # %d, with %d experiments: ', k, N_exp)
    fprintf('\n')
    exp_seed = floor(linspace(1,100,N_exp));            % Specify random seed (determines how the system is excited)
    awgn_seed = round(linspace(1,30,N_exp)); % Random see that determines the form of added noise
    rng(0), snr = randi([35, 70], N_exp,1); % signal to noise ratios
    % rng(0), dt = randi([5, 5],N_exp,1)*0.001;   % Sampling time
    rng(0), dt = randi([5, 10],N_exp,1)*0.001;   % Distinct Sampling time
    fs = 1./dt;                      % Sampling Freq
    rng(0), tend = randi([50, 300], N_exp, 1);   % Simulation Time
    rng(0), Amp = randi([5, 10], N_exp, 1);        % Forcing amplitude
    Y_cell = cell(N_exp,1);          % Create a cell array of candicate experiments
    dY_cell = cell(N_exp,1);         % Create a cell array for dY for all experiments
    % Exogenous Loading function
    close all
    for j = 1:N_exp
        tspan = 0:dt(j):tend(j);           % Discrete time vector
        rr = 1;                           % Define the resampling rate (default is 1)
        rng(exp_seed(j));                    % Specify the seed for consistent results
        loading = 'Random';               % Specify the type of loads
        timing.Tstart = 0;                % Specify the starting time of the simulation
        timing.Tend = tend(j);               % Specify the termination time of the simulation
        timing.dT = dt(j)*rr;                % Samping period (after resampling)
        finfo.A = Amp(j);                     % Excitation amplitude
        finfo.omega = [];                 % Excitation frequency
        finfo.windtype = 'Rectangular';   % Window type if windowing of the excitation is involved
        finfo.pwind = 0;                  % Period of windowing (default is 0, no windowing)
        finfo.floc = 'beam';              % Excitation location (not in use in this example)
        Tend = timing.Tend;
        RandomForcing = [];
        for i = 1:n
            [forcing, ~, tspan, fs(j), dt(j), ~] = forcingfnc(loading, timing, finfo, rr);
            RandomForcing = [RandomForcing, ppval(forcing, tspan)'];
        end

        % Forced oscillation
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
        % Assign the data to the cell arrays
        Y_cell{j} = Y;
        dY_cell{j} = dY;
        % Add Noise to the state variables
        Y = awgn(Y,snr(j),'measured',awgn_seed(j));
        dY = awgn(dY,snr(j),'measured',awgn_seed(j));
    end

    % ------------------------------------------------------------------------
    % Conduct the considered decompositions
    warning('off','all')
    % SAVE SYSTEM INFORMATION AND ============================================
    SaveData = 0; % Data not saved as this stage of dev.
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

    % -------------------------- DECOMPOSITION BEGINS ------------------------
    close all
    mode_indx = [1,4,6,9];
    Metric = 'rell2';
    % ------------------------------------------------------------------------
    % KCVD - GCVD with pre-determined kernal (Left EVP formulation) ==========
    % Fuse multiple data sets into a designed kernel
    C = zeros(2*n, 2*n);
    Y = [];
    dY = [];
    for j = 1:N_exp
        Y_temp = Y_cell{j};
        dY_temp = dY_cell{j};
        K11 = Y_temp'*Y_temp;
        K21 = dY_temp'*Y_temp;
        C = C + 1/2*(K11/K21 + K21/K11);
        Y = [Y; Y_temp];
        dY = [dY; dY_temp];
    end
    C = C/N_exp;

    % Conduct KCVD (left EVP formulation in the paper)
    tic
    [U, V, S, Lambda, W_gcvd] = kcvd(Y,dY,C);
    T_gcvd = toc;
    % Sort the modes according to the frequencies
    [Phi_gcvds, Lambda_gcvds] = ModeSortFreq(inv(W_gcvd'), S{2}/S{1}, 'Uniquedisp');
    % Plot preparation -> decomposing complex modes into real and imag parts
    Phi_gcvd = ReImModes(Phi_gcvds, 'mdof_ff');
    % Visualize the MAC values
    [mac_gcvd,~,~] = MAC(Phi_gcvds(1:10, :), modes(2:end-1, :));
    % Calculate the identification error
    ephi_gcvd = ComplexModeError(modes, [zeros(1,n); Phi_gcvds; zeros(1,n)], Metric);

    % Compare the poles
    mps_gcvd = poles2mp(diag(Lambda_gcvds));
    [~, Lambda_gcvds_full] = ModeSortFreq(inv(W_gcvd'), S{2}/S{1}, 'Full');

    % % Adjust the estimated system poles
    % P1 = U*S{1};
    % % Visualize the adjusted pole estimates
    % Lambda_gcvd_bar = pinv(P1)*((dY-(F*B'))*W_gcvd);
    % [~, Lambda_gcvd_bar_full] = ModeSortFreq(inv(W_gcvd'), Lambda_gcvd_bar, 'Full');
    % xlim([-0.7 0])

    % Caculate the rectified poles
    % mps_gcvd_r = poles2mp(Lambda_gcvd_bar_full);

    % print('-dpng','-r600',strjoin({FigurePath,['Poles_',num2str(n),'dof_complex_mode_comp_gcvd_A_Q3_Q4.png']},filesep))
    % print('-dpng','-r600',strjoin({FigurePath,['Cmodes_',num2str(n),'dof_complex_mode_comp_gcvd_A_Q3_Q4.png']},filesep))

    % ERA ====================================================================
    % ERA will not work if the sampling time among data sets are different

    % ------------------------------------------------------------------------
    % ITD ====================================================================
    % Conduct kernel ITD
    % Compute the kernel for ITD
    % Allocate memory for the Teopeliz matrices
    tic
    Nr = size(Y_cell,1);
    Ncol = size(Y_cell{1},2);
    T11 = zeros(2*Ncol, 2*Ncol);
    T12 = zeros(2*Ncol, 2*Ncol);
    T21 = zeros(2*Ncol, 2*Ncol);
    T22 = zeros(2*Ncol, 2*Ncol);
    ndelay = 3;
    rhank = (ndelay - 1)*Ncol;
    C = zeros(2*Ncol, 2*Ncol);
    % Loop over the number of trials/experiments if mulitple data sets are
    % provided
    for i = 1:Nr
        Y0 = Y_cell{i}; % Extract the ith recording
        if size(Y0,1) > size(Y0,2)
            Y0 = Y0';
        end
        % Hankel Structure Construction (Delay coordinate embedding)
        H = GenHankel(Y0, ndelay);
        H1 = H(1:rhank,:);
        H2 = H(rhank + 1:rhank*2,:);
        % Toeplitz Structure Construction
        T11 = H1*H1' + T11; % T11
        T12 = H1*H2' + T12; % T12
        T21 = H2*H1' + T21; % T21
        T22 = H2*H2' + T22; % T22
    end
    A1 = T21*pinv(T11);
    A2 = T22*pinv(T12);
    C = (A1 + A2)/2;

    % Conduct ITD
    [Lambda_itd, mode_itd] = kITD(C);
    T_itd = toc;
    % Get the ITD continuous-time eigenvalue matrix
    Lambda_itd = diag(log(diag(Lambda_itd))/dt(1)); % Note that the ITD needs the dt to recover the continuous time pole
    % Sort the modes and obtain only the unique displacement modes
    [Phi_itds, Lambda_itds] = modeSortMAC(mode_itd(1:20,:), Lambda_itd, modes(2:end-1,:),'Uniquedisp');

    % Isolate the displacement modes
    [mac_itd,~,~] = MAC(Phi_itds, modes(2:end-1, :));
    Phi_itd = ReImModes(Phi_itds, 'mdof_ff');
    % Calculate the identification error
    ephi_itd = ComplexModeError(modes, [zeros(1,n); Phi_itds(:,1:n); zeros(1,n)], Metric);

    % Compare the poles (this shows that the poles identified by ITD are not necessarily the physical poles of the system)
    mps_itd = poles2mp(diag(Lambda_itds));
    mps_itd = mps_itd(:,1:n);
    [~, Lambda_itds_full] = modeSortMAC(mode_itd(1:20,:), Lambda_itd, modes(2:end-1,:), 'Full');

    % Update metric cells
    % Mode error cell
    ephi_cell_cvd{k} = ephi_gcvd;
    ephi_cell_itd{k} = ephi_itd;
    % MAC cell
    mac_cell_cvd{k} = mac_gcvd;
    mac_cell_itd{k} = mac_itd;
    % MPS cell
    mps_cell_cvd{k} = mps_gcvd;
    mps_cell_itd{k} = mps_itd;

end % End N_exp for loop
disp('All ensemble experiments are done!')
fprintf('\n')
%% Process the mps error data
disp('------ Postprocessing modal parameters error data ------')
[m_, n_] = size(mps);
mps_mat_cvd = zeros(m_, n_, max_N_exp);
mps_mat_itd = zeros(m_, n_, max_N_exp);

for i = 1:max_N_exp
mps_mat_cvd(:,:,i) = mps_cell_cvd{i};
mps_mat_itd(:,:,i) = mps_cell_itd{i};
end
mps_true = repmat(mps, [1,1,max_N_exp]);
mps_err_mat_cvd = abs(mps_mat_cvd - mps_true)./mps_true;
mps_err_mat_itd = abs(mps_mat_itd - mps_true)./mps_true;
%% Visualize the statistics of the natural frequency estimates
% for k = 1:max_N_exp
% figure(k),clf
% plotFreqError({mps_cell_cvd{k}, mps_cell_itd{k}}, mps,...
%     {'o-','^-','v-','x-','*-'}, ...             % Define the markers
%     {}, ...                                % Use default colors
%     [0.2 0.3 0.4 0.5 0.6], ...          % Define alphas
%     {'GCVD','ITD'})      % Define Legends
% set(gca,'yscale','log')
% end
%% Visualize the mode shape errors
% Visualize the damping ratios errors as two subplots as the number of
% trials increases
disp('---------- Visualizing mode shapes errors (fig.1) --------------')
figure(1),clf
DefaultFontSize = 10;
Colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
markers = {'o','^','v','x','*'};
alphas = [0.3 0.4 0.5 0.6 1];
for k = 1:1:max_N_exp
    subplot(1,2,1)
    ephi_cvd = ephi_cell_cvd{k};
    if k == max_N_exp
        p1 = plot(ephi_cvd, 'Color', Colors{1}, 'LineWidth', 1.5);
    else
        p1 = plot(ephi_cvd, 'Color', Colors{1});
    end
    hold on
    scatter(1:10, ephi_cvd, 12, [0 0.4470 0.7410], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1)
    p1.Color =  [p1.Color alphas(k)];
    set(gca, 'yscale', 'log')
%     ylim([5e-5 2e-1])
    xlim([0.8, 10.2])
    xticks(1:10)
    subplot(1,2,2)
    ephi_itd = ephi_cell_itd{k};
    if k == max_N_exp
            p2(k) = plot(ephi_itd, 'Color', Colors{2}, 'linewidth', 1.5);
    else
    p2(k) = plot(ephi_itd, 'Color', Colors{2});
    end
    hold on
    p2(k).Color = [p2(k).Color alphas(k)];
    s2(k)=scatter(1:10, ephi_itd, 12, [0.8500 0.3250 0.0980], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1);
    xlim([0.8, 10.2])
    xticks(1:10)
%     ylim([1 2])
    set(gca, 'yscale', 'log')
end
subplot(1,2,1)
ylabel('$e_\phi$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
GCA = gca;
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(a1)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

subplot(1,2,2)
ylabel('$e_\phi$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
set(gcf,'Paperposition',[0 0 8 1.5])
set(gcf, 'Renderer', 'painters')
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(a2)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

leg = legend(s2, string(n_exp), "NumColumns", 5, "Interpreter","latex",'Location','south');
leg.ItemTokenSize = [13,13]; 
%% Visualize the natural frequency errors
% Visualize the natural frequency errors as two subplots as the number of
% trials increases
disp('---------- Visualizing natural frequency errors (fig.2) --------------')
figure(2),clf
DefaultFontSize = 10;
    Colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
markers = {'o','^','v','x','*'};
alphas = [0.2 0.3 0.4 0.5 1];
for k = 1:max_N_exp
    subplot(1,2,1)
    mps_err_cvd = mps_err_mat_cvd(:,:,k);
    if k == max_N_exp    
        p1 = plot(mps_err_cvd(2,:), 'Color', Colors{1}, 'LineWidth',1.5);
    else
        p1 = plot(mps_err_cvd(2,:), 'Color', Colors{1});
    end
    hold on
    scatter(1:10, mps_err_cvd(2,:), 12, [0 0.4470 0.7410], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1)
    p1.Color =  [p1.Color alphas(k)];
    set(gca, 'yscale', 'log')
    xticks(1:10)
    ylim([5e-5 1e-1])
    xlim([0.8, 10.2])
    subplot(1,2,2)
    mps_err_itd = mps_err_mat_itd(:,:,k);
    if k == max_N_exp
            p2(k) = plot(mps_err_itd(2,:), 'Color', Colors{2}, 'LineWidth',1.5);
    else
        p2(k) = plot(mps_err_itd(2,:), 'Color', Colors{2});
    end
    hold on
    p2(k).Color = [p2(k).Color alphas(k)];
    s2(k)=scatter(1:10, mps_err_itd(2,:), 12, [0.8500 0.3250 0.0980], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1);
    ylim([5e-5 1e-1])
    xlim([0.8, 10.2])
        xticks(1:10)
    set(gca, 'yscale', 'log')
end
subplot(1,2,1)
ylabel('$e_f$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
GCA = gca;
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(b1)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')
subplot(1,2,2)
ylabel('$e_f$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(b2)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')   
set(gcf,'Paperposition',[0 0 8 1.5])
%% Visualize the damping ratio errors
% Visualize the damping ratios errors as two subplots as the number of
% trials increases
disp('---------- Visualizing damping ratio errors (fig. 3)--------------')
figure(3),clf
DefaultFontSize = 10;
Colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
markers = {'o','^','v','x','*'};
alphas = [0.3 0.4 0.5 0.6 1];
for k = 1:1:max_N_exp
    subplot(1,2,1)
    mps_err_cvd = mps_err_mat_cvd(:,:,k);
    if k == max_N_exp    
        p1 = plot(mps_err_cvd(3,:), 'Color', Colors{1}, 'LineWidth',1.5);
    else
        p1 = plot(mps_err_cvd(3,:), 'Color', Colors{1});
    end    
    hold on
    scatter(1:10, mps_err_cvd(3,:), 12, [0 0.4470 0.7410], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1)
    p1.Color =  [p1.Color alphas(k)];
%     set(gca, 'yscale', 'log')
    ylim([5e-5 3e-1])
    xlim([0.8, 10.2])
            xticks(1:10)

    subplot(1,2,2)
    mps_err_itd = mps_err_mat_itd(:,:,k);
    if k == max_N_exp    
        p2(k) = plot(mps_err_itd(3,:), 'Color', Colors{2}, 'LineWidth',1.5);
    else
        p2(k) = plot(mps_err_itd(3,:), 'Color', Colors{2});
    end
    hold on
    p2(k).Color = [p2(k).Color alphas(k)];
    s2(k)=scatter(1:10, mps_err_itd(3,:), 12, [0.8500 0.3250 0.0980], ...
        markers{k}, 'alphadata', alphas(k), ...
        'MarkerEdgeAlpha', alphas(k), 'LineWidth', 1);
    xlim([0.8, 10.2])
    ylim([1.3 1.8])
            xticks(1:10)

%     set(gca, 'yscale', 'log')
end
subplot(1,2,1)
ylabel('$e_\zeta$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
GCA = gca;
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(c1)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

subplot(1,2,2)
ylabel('$e_\zeta$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
set(gcf,'Paperposition',[0 0 8 1.5])
set(gcf, 'Renderer', 'painters')
text(GCA.Position(1)-0.3,GCA.Position(2) + GCA.Position(4)+0.05, '(c2)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')


%% Visualize error statistics using bar plots
disp('---------- Visualizing error plots using bar plots (fig.4) ---------------')
figure(4),clf
subplot(1,3,1)
bar([mean(ephi_cell_cvd{1}),mean(ephi_cell_itd{1});mean(ephi_cell_cvd{2}),mean(ephi_cell_itd{2});mean(ephi_cell_cvd{3}),mean(ephi_cell_itd{3});mean(ephi_cell_cvd{4}),mean(ephi_cell_itd{4});mean(ephi_cell_cvd{5}),mean(ephi_cell_itd{5})])
% xlabel('Number of total experiments', 'interpreter','latex')
ylabel('$e_\phi$','Interpreter','latex')
xticks(1:max_N_exp)
xticklabels(string(n_exp))
set(gca, 'yscale', 'log')
ylim([1e-3 0.5])
yticks([1e-3, 1e-2, 1e-1])
leg = legend({'GCVD','ITD'},'NumColumns',2,'location','northwest','interpreter','latex');
leg.ItemTokenSize = [13,13]; 
GCA = gca;
text(GCA.Position(1)-0.35,GCA.Position(2) + GCA.Position(4)+0.05, '(d1)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

subplot(1,3,2)
bar([mean(mps_err_mat_cvd(2,:,1)),mean(mps_err_mat_itd(2,:,1));mean(mps_err_mat_cvd(2,:,2)),mean(mps_err_mat_itd(2,:,2));mean(mps_err_mat_cvd(2,:,3)),mean(mps_err_mat_itd(2,:,3));mean(mps_err_mat_cvd(2,:,4)),mean(mps_err_mat_itd(2,:,4));mean(mps_err_mat_cvd(2,:,5)),mean(mps_err_mat_itd(2,:,5))])
set(gca, 'yscale', 'log')
xlabel('Number of total experiments','Interpreter','latex')
ylabel('$e_f$','Interpreter','latex')
xticks(1:max_N_exp)
xticklabels(string(n_exp))
yticks([1e-3 1e-2 1e-1])
yticklabels({'$10^{-3}$','$10^{-2}$','$10^{-1}$'})
ylim([1e-3 0.2])
text(GCA.Position(1)-0.4,GCA.Position(2) + GCA.Position(4)+0.05, '(d2)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

subplot(1,3,3)
bar([mean(mps_err_mat_cvd(3,:,1)),mean(mps_err_mat_itd(3,:,1));mean(mps_err_mat_cvd(3,:,2)),mean(mps_err_mat_itd(3,:,2));mean(mps_err_mat_cvd(3,:,3)),mean(mps_err_mat_itd(3,:,3));mean(mps_err_mat_cvd(3,:,4)),mean(mps_err_mat_itd(3,:,4));mean(mps_err_mat_cvd(3,:,5)),mean(mps_err_mat_itd(3,:,5))])
ylim([0.02 2])
set(gca, 'yscale', 'log')
% xlabel('Number of total experiments','Interpreter','latex')
ylabel('$e_\zeta$','Interpreter','latex')
yticks([0.01, 0.1, 1])
yticklabels({'$10^{-2}$','$10^{-1}$','$10^0$'})
xticks(1:max_N_exp)
xticklabels(string(n_exp))
set(gcf,'paperposition',[0,0,8,1.5])
text(GCA.Position(1)-0.35,GCA.Position(2) + GCA.Position(4)+0.05, '(d3)', 'Units','normalized','FontSize',DefaultFontSize, 'Interpreter','latex')

disp('Emsemble Mode-wise Averaged Mode Error:')
disp('    GCVD      ITD')
disp([mean(ephi_cell_cvd{1}),mean(ephi_cell_itd{1});mean(ephi_cell_cvd{2}),mean(ephi_cell_itd{2});mean(ephi_cell_cvd{3}),mean(ephi_cell_itd{3});mean(ephi_cell_cvd{4}),mean(ephi_cell_itd{4});mean(ephi_cell_cvd{5}),mean(ephi_cell_itd{5})])
disp('Emsemble Mode-wise Averaged Natural Frequency Error:')
disp('    GCVD      ITD')
disp([mean(mps_err_mat_cvd(2,:,1)),mean(mps_err_mat_itd(2,:,1));mean(mps_err_mat_cvd(2,:,2)),mean(mps_err_mat_itd(2,:,2));mean(mps_err_mat_cvd(2,:,3)),mean(mps_err_mat_itd(2,:,3));mean(mps_err_mat_cvd(2,:,4)),mean(mps_err_mat_itd(2,:,4));mean(mps_err_mat_cvd(2,:,5)),mean(mps_err_mat_itd(2,:,5))])
disp('Emsemble Mode-wise Averaged Damping Ratio Frequency Error:')
disp('    GCVD      ITD')
disp([mean(mps_err_mat_cvd(3,:,1)),mean(mps_err_mat_itd(3,:,1));mean(mps_err_mat_cvd(3,:,2)),mean(mps_err_mat_itd(3,:,2));mean(mps_err_mat_cvd(3,:,3)),mean(mps_err_mat_itd(3,:,3));mean(mps_err_mat_cvd(3,:,4)),mean(mps_err_mat_itd(3,:,4));mean(mps_err_mat_cvd(3,:,5)),mean(mps_err_mat_itd(3,:,5))])
%%
disp('===================================================================')