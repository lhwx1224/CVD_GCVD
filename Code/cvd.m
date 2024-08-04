% Characteristic Value Decomposition
% Conventions: for the input spatiotemporal field measurement matrix X, time
% evolves row-wise, spatial variable populates column-wise. If the matrix
% size is m-by-n, then, there are m time samples and n spatial samples.

function [U, Sigma, V, W] = cvd(X, dX, Method, dt)
if nargin == 1
    dX = diff(X);
end

if nargin < 3
    Method = 'projective';
end

[m, n] = size(X);            % Get the dimensions of the input data
% Caclulate the cross-covariance matrices
Kxx = X'*X;                  % Covariance of the field measurement matrix
Kxdx = X'*dX;                % Cross-covariance between the field measurement and the derivative
% -------------- Exploratory Code ------------------
Kdxx = dX'*X;             % Debug purpose
% --------------------------------------------------
% Conduct the generalized eigenvalue problem
% [V, Sigma, W] = eig(Kxdx, Kxx);  % Conduct generalized eigenvalue decomposition
% Left eigenvector formulation
[V, Sigma, W] = eig(Kdxx, Kxx);  % Conduct generalized eigenvalue decomposition
% Use the following eigenvalue problem to yield cvd mode as the right
% eigenvectors
% [V, Sigma, W] = eig(Kdxx/Kxx);  % Conduct generalized eigenvalue decomposition

% Calculate the temporal coordinates
if isequal(lower(Method), 'projective')   % IF projective method
    U = X*W;
elseif isequal(lower(Method), 'constructive')  % IF constructive method
    t = 0:dt:(m-1)*dt;
    T = repmat(t,n,1);
    eta0 = V*X(1,:)'; % Modal coordinate 
    U = exp(Sigma*T)*diag(eta0); % Constructed modal coordinates
    U = U';
end
end

% t = 0:0.01:100;
% tt = repmat(t, 2, 1);
% Lambda = [-0.1+3.14i, 0; 0, -0.2+3.14*5i];
% plot(imag(exp(Lambda*tt)'*[1 0; 0 2]))