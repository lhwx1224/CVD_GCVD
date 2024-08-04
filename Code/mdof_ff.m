function [M, C, K] = mdof_ff(n, m, c, k, index_Ic, ci)
% A one-degree-of-freedom oscillator chain w/ both ends fixed
% input: n - number of degree of freedoms
%        m - mass vector dictates the masses of the point masses
%        c - damping vector dictates the damping coeffs. of the n+1 dampers
%        k - stiffness vector dicates the stiffness of the n+1 springs
%        index_Ic - index vector of the additional dampers
%        ci - (sparse) damping vector of the additional dampers
% 
% output: M - mass matrix
%         C - damping matrix
%         K - stiffness matrix
%
% Hewenxuan Li 2023-04-08 @ Cornell

% System Parameters
% check inputs
if nargin < 2
    m = ones(n,1);            % mass vector (use 1 if not provided)
    k = 4*ones(n+1,1);        % stiffness vector (use 1 if not provided)
    c = 0.3*ones(n+1,1);      % damping vector (use 0.3 if not provided)
    index_Ic = [1 2];         % index of the additional dampers
    ci = 0.3*ones(n+1,1);     % viscous damping coeffs. of the additional dampers
elseif nargin < 5
    index_Ic = [1 2];         % index of the additional dampers
    ci = 0.3*ones(n+1,1);     % viscous damping coeffs. of the additional dampers
end

Ic = zeros(n+1,1);        % initialize indicator to the additinal dampers (for non-proportional damping)
Ic(index_Ic) = 1;         % input the index of the additional dampers

% Assemble the parameters into system matrices
M = zeros(n);             % initialize mass matrix
C = zeros(n);             % initialize damping matrix
K = zeros(n);             % initialize stiffness matrix

%% Left boudary point mass
for i = 1
    % Mass matrix
    M(i,i) = m(i);
    % Damping matrix
    C(i,i) = c(i) + c(i+1);
    C(i, i+1) = -c(i+1);
    % Additional dampers
    if Ic(i) == 1
        C(i, i) = C(i,i) + ci(i);
    end
    % Stiffness matrix
    K(i,i) = k(i) + k(i+1);
    K(i,i+1) = -k(i+1);
end

%% Right boudary point mass
for i = n
    % Mass matrix
    M(i,i) = m(i);
    % Damping matrix
    C(i,i) = c(i) + c(i+1);
    C(i, i-1) = -c(i);
    % Additional dampers
    if Ic(i) == 1
        C(i, i) = C(i,i) + ci(i);
        C(i, i-1) = C(i,i) - ci(i);
    end
    % Stiffness matrix
    K(i,i) = k(i) + k(i+1);
    K(i,i-1) = -k(i);    
end

%% Inner point masses
if n > 2
    for i = 2:n-1
        % Mass matrix
        M(i,i) = m(i);
        % Damping matrix
        C(i,i) = c(i) + c(i+1);
        C(i, i+1) = -c(i+1);
        C(i, i-1) = -c(i);
        % Additional dampers
        if Ic(i) == 1
            C(i, i) = C(i,i) + ci(i);
            C(i, i-1) = C(i,i) - ci(i);
        end
        % Stiffness matrix
        K(i,i) = k(i) + k(i+1);
        K(i,i+1) = -k(i+1);
        K(i,i-1) = -k(i);
    end
end