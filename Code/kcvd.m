function [U, V, S, Lambda, W, R] = kcvd(X1, X2, C)
% KCVD.m conducts the *kernel characteristic value decomposition* to the
% given set of two matrices whose shape are identical. 
% Syntax: [U, V, S, Lambda, W, R] = gcvd(X1,X2)
% 
% Input:  C - Kernel matrix to be decomposed via CVD
% Output: U - Common left basis matrix (cooridnates with unit column norm)       
%         V - Common right basis matrix (modes)
%         S - a cell array with real and complex magnitude matrices for the
%             first and the second data matrices respectively
%         Lambda - Complex eigenvalues to the balance eigenvalue problem
%         R - complex characteristic value matrix to the GCVD problem 
%
% Created by Hewenxuan Li, June 19, 2021
% Modifications            August, 2023 - forced case option added
%                          September, 2023 - changed the notations, updated
%                          the code descriptions
% -------------------------------------------------------------------------

% Check the data size and rank
[m, n] = size(C);
if m ~= n
    error('Input kernel matrix C has to be square!')
end

% CVD
[~, Lambda, W] = eig(C);    % Eigendecomposition to the balance quotient matrix
V = inv(W');                   % Calculate the mode matrix (adjoint of the projective modes)
P1 = X1*W;                     % Projected data matrix (projective coordinates)
P2 = X2*W;                     % Projected velocity matrix (projective rate)
S1 = diag(sqrt(sum(abs(P1).^2))); % Magnitude matrix to Y1
U = P1/S1;                 % Common unit-norm coordinate matrix
S2 = S1*pinv(P1)*P2;       % Complex magnitude to the second data matrix via U = M1*inv(Sigma1) = M2*inv(Sigma2) -> Sigma2 = Sigma1*pinv(M1)*M2
S = {S1, S2};              % Assign the magntiudes into S
R = S{2}/S{1};             % Characteristic value matrix
end
