function Err = ComplexModeError(Phi, Phi_hat, Metric)
% ComplexModeError(Phi, Phi_hat, Metric) Calculate the error between the
% estimated mode and the true mode.
% 
% syntax: Err = ComplexModeError(Phi, Phi_hat, Metric)
% input:  Phi - 2d double, the true mode shapes
%         Phi_hat - 2d double, the estimated mode shapes
%         Metric - string, the type of metrics used to characterize the
%         error
% output: Err - double, the error metric calculated
%
% Created by Hewenxuan Li @ Cornell 2023-04-04
% Modified on 2023-04-24, changed the mode error metrics all depend on
% complex-valued input modes.
% Added the subspace angle option - Method = 'angle'
% Added the rmse option - Method = 'rmse'

% Check if the modes are complex-valued
n = size(Phi, 2);
if nargin < 3
    Metric = 'rell2';
end
Normalization = 'norm';
% Normalize the complex mode shapes
Phi_n = normalize(Phi, Normalization);
Phi_hat_n = normalize(Phi_hat, Normalization);

if isreal(Phi)
    warning('Input modes should be complex valued, ModeError() will be called!')
    Err = ModeError(Phi, Phi_hat, Metric);
else
    Err = zeros(n,1); % Allocate memory for the error metric
    if isequal(lower(Metric),'rell2') % Relative L2 error
        % Phase matching
        Phi_hat_n_pm = ModePhaseMatch(Phi_hat_n, Phi_n, n, 'complex');
        % Check if all modes are properly matched
        indx = find(sum(isnan(Phi_hat_n_pm),1),n+2);
        if ~isempty(indx)
            warning(['Phase-mathed modes failed! mode index:', num2str(indx), ' are detected to be NaNs!'])
            warning('NaN valued modes are replaced by the original mode estimates. Please consider using other metrics!')
            Phi_hat_n_pm(:,indx) = Phi_hat_n(:,indx);
        end
        % Calculate the error
        dPhi = Phi_hat_n_pm - Phi_n;
        for i = 1:n
            Err(i) = (dPhi(:,i)'*dPhi(:,i))/(Phi_n(:,i)'*Phi_n(:,i));
        end
    elseif isequal(lower(Metric),'angle') % Subspace angle
        % Subspace angle
        for i = 1:n
            Err(i) = subspace(Phi_hat_n(:,i), Phi_n(:,i));
        end
    elseif isequal(lower(Metric),'rmse') % 
        % Phase matching
        Phi_hat_n_pm = ModePhaseMatch(Phi_hat_n, Phi_n, n, 'complex');
        % Calculate the error
        dPhi = Phi_hat_n_pm - Phi_n;
        % For each mode shape
        for i = 1:n
            Err(i) = sqrt(mean(abs(dPhi(:,i)).^2)); % Calculate the RMSE
        end
    end
end

