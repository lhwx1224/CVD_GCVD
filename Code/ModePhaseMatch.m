function Phi_hat_pm = ModePhaseMatch(Phi_hat, Phi, n, method)
% ModePhaseMatch matches the phases between the estimated modes 'Phi_hat'
% and the true modes 'Phi' for 'n' number of modes; the phase-matched mode
% shapes are denoted as 'Phi_hat_pm'. The method is based on a
% psudo-inverse as a rotational transformation to the estimated mode shape.
% syntax:  Phi_hat_pm = ModePhaseMatch(Phi_hat, Phi, n)
% input: Phi_hat - Estimated mode shapes
%        Phi - true mode shapes
%        n - number of modes to match
% ouput: Phi_hat_pm - phase-matched estimated mode shapes
% Hewenxuan Li 04-04-2023 @ Cornell

% Define the phase-matched mode vectors
if nargin < 3 
    n = size(Phi,1);
    method = 'ReIm';
elseif nargin < 4
    method = 'ReIm';
end

if isequal(lower(method), 'reim')
Phi_hat_pm = cell(n,1);
for k = 1:n
    PT = ([Phi_hat{k}(:,1) Phi_hat{k}(:,2)]'*[Phi_hat{k}(:,1) Phi_hat{k}(:,2)])\...
         [Phi_hat{k}(:,1) Phi_hat{k}(:,2)]'*[real(Phi(:,k)) imag(Phi(:,k))]; % Phase matching transform, psudo-inverse based
    rim = [Phi_hat{k}(:,1) Phi_hat{k}(:,2)]*PT;                              % 
    Phi_hat_pm{k} = [rim(:,1), rim(:,2)];
    if sum(isnan(Phi_hat_pm{k}),'all') ~= 0
        Phi_hat_pm{k} = Phi_hat{k};
    end
end
elseif isequal(lower(method), 'complex')
Phi_hat_pm = zeros(size(Phi));
for k = 1:n
    PT = ([real(Phi_hat(:,k)), imag(Phi_hat(:,k))]'*[real(Phi_hat(:,k)), imag(Phi_hat(:,k))])\...
         [real(Phi_hat(:,k)), imag(Phi_hat(:,k))]'*[real(Phi(:,k)) imag(Phi(:,k))]; % Phase matching transform, psudo-inverse based
    rim = [real(Phi_hat(:,k)) imag(Phi_hat(:,k))]*PT;                              % 
    Phi_hat_pm(:,k) = rim(:,1) + 1i*rim(:,2);
%     if sum(isnan(Phi_hat_pm{k}),'all') ~= 0
%         Phi_hat_pm{k} = Phi_hat{k};
%     end
end
end
