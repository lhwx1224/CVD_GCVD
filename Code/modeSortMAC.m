function [modes_s, poles_s] = modeSortMAC(modes_hat, poles_hat, modes, varargin)

% Check the output type (default, JordanPair)
if nargin > 3
    type = varargin{1};
end
% check if the poles are vectors (if a matrix, reshape it into a vector)
if size(poles_hat, 2) ~= 1 || size(poles_hat, 1) ~= 1
    poles_hat = diag(poles_hat);
end
[Nr_, Nc_] = size(modes);            % Nr_, Nc_ number of rows and columns of the true modes
[Nr, Nc] = size(modes_hat);              % Nc - number of estimated modes (including complex conjugates)
Nm = floor(Nc/2);                        % Nm = Nc/2 - number of possible complex modes;
Phi = zeros(Nr/2,Nm);                    % mode shape matrix (different from the sort by frequency code)

% Make the poles to a designated order with positive imag part appear
% first
for i = 1:Nm
    if imag(poles_hat(2*i-1)) < imag(poles_hat(2*i))
        temp_pole = poles_hat(2*i-1);
        poles_hat(2*i-1) = poles_hat(2*i);
        poles_hat(2*i) = temp_pole;
    end
end

lambda = poles_hat(1:2:end);             % Extract the unique poles

for i = 1:Nm
    ii = 2*i - 1;
    Phi(:,i) = modes_hat(1:Nr/2,ii); % Mode shape
end

% Conduct MAC and output the sort index
[~, ~, Is] = MAC(Phi, modes); % Conduct MAC between the true modes and the identified unique modes

if nargin == 2 % Default output of modes
    modes_s = Phi(:,Is);                           % Nc*Nc/2 sorted modes (no conjugate pair)
    poles_s = [lambda(Is); conj(lambda(Is))];       % Return poles in two columns of complex eigenvalues
elseif isequal(lower(type), 'unique')
    modes_s = Phi(1:Nc,Is);
    poles_s = lambda(Is);                          % Return the unique poles (upper half in the complex plane)
elseif isequal(lower(type), 'uniquedisp')
    % Return only the Nm*Nm complex mode shapes which correspond to the
    % modes of the displacement response
    modes_s = Phi(1:Nr_,Is);
    poles_s = lambda(Is);                          % Return the unique poles (upper half in the complex plane)
elseif isequal(lower(type), 'full')
    modes_s = zeros(size(modes_hat));
    poles_Nm = lambda(Is); % sort Nm-number unique poles
    poles_s = zeros(size(poles_hat)); % Initiate full size sorted poles
    % Assign the full sorted modes and full sorted poles
    for i = 1:length(Is)
        modes_s(:,i*2-1) = modes_hat(:,Is(i)*2-1); % Complex modes sorted
        modes_s(:,i*2) = conj(modes_hat(:,Is(i)*2-1)); % Conjugate modes sorted
        poles_s(2*i-1,2*i-1) = poles_Nm(i); % Comples poles sorted
        poles_s(2*i,2*i) = conj(poles_Nm(i)); % Conjugate poles sorted
    end
end

% TODO:
% The following code is from the ModeSortFreq, which can be incorporated
% into a modeSort function in the future
% 
% % Check if the modal paramters are complex
% if ~isreal(poles_hat)
%     [Nr, Nc] = size(modes_hat);                  % Nm - number of modes
%     Nm = floor(Nc/2);                        % Nm = Nc/2 - number of complex modes;
%     omega_hat = zeros(1,Nm);                 % Omegan - natural frequencies
%     zeta_hat = zeros(1,Nm);                  % Zeta - modal damping freqs
%     Phi = zeros(Nr,Nm);                      % mode shape matrix
%     for i = 1:Nm
%         ii = 2*i - 1;    % Get one of the characteristic value
%         omega_hat(i) = abs(poles(ii)); % Absolute value to the characteristic value
%         zeta_hat(i) = - real(poles(ii))/omega_hat(i); % -real(lambda)/abs(lambda) (Normalized real part)
%         Phi(:,i) = modes(:,ii); % Mode shape
%     end
%     [omega_hat_sort, I] = sort(omega_hat, 'ascend');
%     zeta_hat_sort = zeta_hat(I);
%     
%     if nargin == 2 % Default output of modes
%         modes_s = Phi(:,I);                  % Nc*Nc/2 sorted modes (no conjugate pair)
%         mp = [omega_hat_sort; zeta_hat_sort]; % modal parameter matrix
%         % Return poles in two columns of complex eigenvalues
%         poles_s = mp2poles(mp);              
%     elseif isequal(lower(type), 'jordanpair')
%         % Nc*2Nc complex mode sorted according to the oscillation freq
%         % First Nc*Nc are the modes related to displacement
%         % Second Nc*Nc are the conjugate modes
%         modes_s = [Phi(1:Nc/2,I), conj(Phi(1:Nc/2,I))];
%         mp = [omega_hat_sort; zeta_hat_sort];
%         poles_s = mp2poles(mp);
%         % 2Nc*Nc complex eigenvalues
%         J = [diag(poles_s(:,1)), zeros(Nm,Nm);...
%             zeros(Nm,Nm), diag(poles_s(:,2))];
%         % Reassign the poles with the eigenvalue matrix J
%         poles_s = J;
%     elseif isequal(lower(type), 'unique')
%         % Return only the Nm*Nm complex mode shapes which correspond to the
%         % modes of the displacement response
%         modes_s = Phi(1:Nc,I);
%         % Return the unique complex eigenvalues by rejecting the conjugate
%         % pairs, Nm*1
%         mp = [omega_hat_sort; zeta_hat_sort];
%         poles_s = mp2poles(mp);
%         poles_s = poles_s(:,1);
%     elseif isequal(lower(type), 'uniquedisp')
%         % Return only the Nm*Nm complex mode shapes which correspond to the
%         % modes of the displacement response
%         modes_s = Phi(1:Nc/2,I);
%         % Return the unique complex eigenvalues by rejecting the conjugate
%         % pairs, Nm*1
%         mp = [omega_hat_sort; zeta_hat_sort];
%         poles_s = mp2poles(mp);
%         poles_s = poles_s(:,1);
%     elseif isequal(lower(type), 'full')
%         modes_s_full = Phi(:,I);
%         modes_s = zeros(size(modes));
%         
%         poles_Nm = zeros(Nm,1); % Initiate Nm-number unique poles
%         % Obtain the complex poles with positive imaginary part
%         for i = 1:Nm
%             poles_Nm(i) = poles(2*i-1); % Assign the unique poles
%         end
%         poles_Nm = poles_Nm(I); % Sort the unique poles
%         poles_s = zeros(size(poles)); % Initiate full size sorted poles
%         % Assign the full sorted modes and full sorted poles
%         for i = 1:Nm
%             modes_s(:,i*2-1) = modes_s_full(:,i); % Complex modes sorted
%             modes_s(:,i*2) = conj(modes_s_full(:,i)); % Conjugate modes sorted
%             poles_s(2*i-1,2*i-1) = poles_Nm(i); % Comples poles sorted 
%             poles_s(2*i,2*i) = conj(poles_Nm(i)); % Conjugate poles sorted
%         end
%     end
% end
% 
