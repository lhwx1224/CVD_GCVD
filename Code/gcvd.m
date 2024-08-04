function [U, V, S, Lambda, W, R] = gcvd(X1,X2,options)
% GCVD.m conducts the *generalized characteristic value decomposition* to the
% given set of two matrices whose shape are identical. The decomposition is
% of the following form:
%                       Y1 = U*S1*V^H
%                       Y2 = U*S2*V^H
% where U and V are the common left basis functions (coordinates) and the
% right basis functions (modes), respectively. U and V are generally
% complex valued matrices. Sigma1 is generically a real valued diagonal
% matrix and S2 is of complex values; the fraction of their diagonal
% element reveals the characteristic values of the underlying dynamics of
% the data (and/or the system which generates the data).
% This algorithm solves for the common mode and common coordinate through a
% balanced eigenvalue decomposition of the sum of the quotients of
% covariance matrices (both cross-covariances and the auto-covariances) of
% the following form:
%
%                       W^H*S = Lambda*W^H
% where
%                       S = 1/4*(K11 + K12 + K21 + K22) -> 
%                           sum of quotients of covariance matrices.
%                       K11 = X1'*X1 -> auto-covariance of X1
%                       K12 = X1'*X2 -> cross-covariance of X1 and X2
%                       K21 = X2'*X1 -> cross-covariance of X2 and X1
%                       K22 = X2'*X2 -> auto-covariance of X2
%
% Syntax: [U, V, S, Lambda, W, R] = gcvd(X1,X2)
% 
% Input:  X1 - first data matrix 
%         X2 - second data matrix
%         options - options to choose which form of the quotient to use
%                   - 'free': all four quotients included
%                   - 'forced': two quotients with pseudoinverses included
%                   - 'forward': two quoteints correspond to forward
%                   dynamics included
%                   - 'backward': two quoteints correspond to backward
%                   dynamics included
%                   - 'cvd': only one quotient that correspond to the CVD
%                   formulation is included
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
if nargin < 3
    options = 'free';
end
% Check the data size and rank
if size(X1,1) ~= size(X2,1) || size(X1,2) ~= size(X2,2)
    [mA, nA] = size(X1);       % Check the size of the first matrix
    [mB, nB] = size(X2);       % Check the size of the second matrix
    mmin = min(mA,mB);         % Check the minimum sizes of the data matrices
    nmin = min(nA,nB);
    r = rank([X1; X2]);        % Rank of the concatenated matrix
    if nmin > r
        X1 = X1(1:mmin,1:r);   % Change the size of the matrix
        X2 = X2(1:mmin,1:r);   
        disp(['Since A and B must have the same shape, they are reshaped to: ',num2str(size(X1,1)),'*',num2str(size(X1,2))])
        disp(['The matrices are rank-deficient, the number of columns has beem adjusted to:',num2str(r)])
    else
        X1 = X1(1:mmin,1:nmin);
        X2 = X2(1:mmin,1:nmin);
        disp(['Since A and B must have the same shape, they are reshaped to: ',num2str(size(X1,1)),'*',num2str(size(X1,2))])
    end
end

% If the rank and size are satisfied or justified
% Construct the cross-covariance matrices
K11 = X1'*X1;                  % K11
K12 = X1'*X2;                  % K12
K21 = X2'*X1;                  % K21
K22 = X2'*X2;                  % K22

% Construct the quotient matrix
if isequal(lower(options), 'free')
    C = 1/4*(K11/K21 + K21/K11 + K12/K22 + K22/K12); % Balanced covariance matrix
elseif isequal(lower(options), 'forced')
    C = 1/2*(K11/K21 + K21/K11); % Balanced quotient matrix wihtout standard pseudoinverses
elseif isequal(lower(options), 'forward')
    C = 1/2*(K21/K11 + K12/K22); % Balanced forward covariance matrix
elseif isequal(lower(options), 'backward')
    C = 1/2*(K11/K21 + K22/K12); % Balanced backward covariance matrix
elseif isequal(lower(options), 'cvd')
    C = K11/K21; % CVD covariance matrix
end

% CVD
[~, Lambda, W] = eig(C);    % Eigendecomposition to the balance quotient matrix
% Modified on Sep 8th, 2023
% W = inv(V');                 % Porjective modes (adjoint modes) to the balanced eigenvalue problem
% P1 = X1/V';                  % Projected data matrix
% P2 = X2/V';                  % Projected velocity matrix
V = inv(W');                   % Calculate the mode matrix (adjoint of the projective modes)
P1 = X1*W;                     % Projected data matrix (projective coordinates)
P2 = X2*W;                     % Projected velocity matrix (projective rate)
S1 = diag(sqrt(sum(abs(P1).^2))); % Magnitude matrix to Y1
U = P1/S1;                 % Common unit-norm coordinate matrix
S2 = S1*pinv(P1)*P2;       % Complex magnitude to the second data matrix via U = M1*inv(Sigma1) = M2*inv(Sigma2) -> Sigma2 = Sigma1*pinv(M1)*M2
S = {S1, S2};              % Assign the magntiudes into S
R = S{2}/S{1};             % Characteristic value matrix
end
