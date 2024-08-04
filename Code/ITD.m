function [poles, modes] = ITD(Y)
% Ibrahim Time Domain Modal Identification
% [poles, modes] = ITD(Y) returns the discrete time poles (natural
% frequencies and dampling ratios) and the modes from a give set of
% responses (assumed to be free decay data sets).
% 
% syntax: [poles, modes] = ITD(Y)
%
% input: Y - cell array, system response matrix (assuming free decays)
%
% output: poles - 2d array, system poles in a matrix form
%         modes - 2d array, system modes in a matrix form
% 
% Number of records (response matrices from sensor output)
% Hewenxuan Li, Output-Only Modal Analysis Toolbox v0.0
% Created 2021 @ URI
% Last modified 09-20-2023 @ Cornell

% Check the data orientation, the default recognizes rows as different
% responses
if size(Y,1) < size(Y,2)
    Y = Y';
end
% Number of recordings
Nr = size(Y,1);
% Size of the recordings (n denotes the size of the state/configuration space)
if size(Y{1},1) > size(Y{1},2)
    [~, n] = size(Y{1});
else
    [n, ~] = size(Y{1});
end

% Allocate memory for the Teopeliz matrices
T11 = zeros(2*n, 2*n);
T12 = zeros(2*n, 2*n);
T21 = zeros(2*n, 2*n);
T22 = zeros(2*n, 2*n);
% Number of delays is set to 3
% ndelay = 2; (legacy)
% rhank = ndelay*n; (legacy)
ndelay = 3;
rhank = (ndelay - 1)*n;
% Loop over the number of trials/experiments if mulitple data sets are
% provided
for i = 1:Nr
    Y0 = Y{i}; % Extract the ith recording
    if size(Y0,1) > size(Y0,2)
        Y0 = Y0';
    end
    if size(Y0,1) ~= n
        error(['The input matrices #',num2str(i),'has incompatible number of sensor inputs!'])
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

A = (A1 + A2)/2;
[EVec, Eval] = eig(A);
poles = sqrt(Eval);
modes = EVec;