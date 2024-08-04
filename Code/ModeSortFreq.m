function [modes_s, poles_s] = ModeSortFreq(modes, poles, varargin)
% MODESORTFREQ sorts the system poles (eigenvalues) and the modes
% (eigenvectors) according to the *DAMPED NATURAL FREQUENCIES* (imaginary
% part of the complex eigenvalues).
%
% Syntax: [modes_s, poles_s] = ModeSortFreq(modes, poles, varargin)
%
% Input: modes - modal matrix obtained from eigen-decomposition
%        poles - eigenvalue matrix obtained from eigen-decomposition
%        varargin - options that designates the type of the returned sorted
%                   modal paramters:
%        if not specified: the returned modes only contains 2*Nm-by-Nm
%        elements, whose columns are the complex mode shape sorted by the
%        angular frequency (imaginary part of the eigenvalues),
%
%              modes_s = [Phi_s; Phi_s*Lambda_s] \in C^{2Nm*Nm}
%              poles_s = [Lambda_s, conj(Lambda_s)] \in C^{Nm*2}
%
%        if varargin{1} == 'JordanPair': returned modes is of Nm*2Nm, whose
%        first Nm columns are the Nm*Nm complex modes and whose second Nm
%        columns are the Nm*Nm complex conjugate pair modes of the first Nm
%        columns,
%
%              modes_s = [Phi_s, conj(Phi_s)] \in C^{Nm*2Nm}
%              poles_s = [Lambda_s, zeros(Nm); 
%                         zeros(Nm), conj(Lambda_s)] \in C^{2Nm*2Nm}
%
%        if varargin{1} == 'Unique': returned modes are of Nm*Nm complex
%        values which are the modes to the underlying dynamic problem. And
%        the returned poles are unique complex eigenvalues of Nm*1 sorted
%        by the angular frequencies.
%
%              modes_s = [Phi_s] \in C^{Nm*Nm}
%              poles_s = [Lambda_s] \in C^{Nm*1}
%
%        if varargin{1} == 'Full': return sorted full-size complex modes of
%        the size 2*Nm-by-2*Nm and their corresponding poles of the same
%        size.
%
% Hewenxuan Li July 14, 2021
% Modified 7/19/2021 Added full-size sorted mode function

% Check the output type (default, JordanPair)
if nargin > 2
    type = varargin{1};
end
% check if the poles are vectors (if a matrix, reshape it into a vector)
if size(poles, 2) ~= 1 || size(poles, 1) ~= 1
    poles = diag(poles);
end

% Check if the modal paramters are complex
if ~isreal(poles)
    [Nr, Nc] = size(modes);                  % Nm - number of modes
    Nm = floor(Nc/2);                        % Nm = Nc/2 - number of complex modes;
    omega_hat = zeros(1,Nm);                 % Omegan - natural frequencies
    zeta_hat = zeros(1,Nm);                  % Zeta - modal damping freqs
    phi = zeros(Nr,Nm);                      % mode shape matrix
    
    % Make the poles to a designated order with positive imag part appear
    % first
    for i = 1:Nm
        if imag(poles(2*i-1)) < imag(poles(2*i))
            temp_pole = poles(2*i-1);
            poles(2*i-1) = poles(2*i);
            poles(2*i) = temp_pole;
        end
    end

    % Sort the modes and poles according to the natural frequencies
    % identified
    if isequal(lower(type), 'full')
        for i = 1:Nc
            omega_hat(i) = abs(poles(i)); % Absolute value to the characteristic value
            zeta_hat(i) = - real(poles(i))/omega_hat(i); % -real(lambda)/abs(lambda) (Normalized real part)
            phi(:,i) = modes(:,i); % Mode shape
        end
        [~, I] = sort(omega_hat, 'ascend');
        poles_s = diag(poles(I));

        % Initialize sorted modes
        modes_s = phi(:,I);
%         modes_s = zeros(size(modes));

%         poles_Nm = zeros(Nm,1); % Initiate Nm-number unique poles
%         % Obtain the complex poles with positive imaginary part
%         for i = 1:Nm
%             poles_Nm(i) = poles(2*i-1); % Assign the unique poles
%         end
%         poles_Nm = poles_Nm(I); % Sort the unique poles
%         poles_s = zeros(Nc,Nc); % Initiate full size sorted poles
%         % Assign the full sorted modes and full sorted poles
%         for i = 1:Nm
%             modes_s(:,i*2-1) = modes_s_full(:,i); % Complex modes sorted
%             modes_s(:,i*2) = conj(modes_s_full(:,i)); % Conjugate modes sorted
%         end
    else
        for i = 1:Nm
            ii = 2*i - 1;    % Get one of the characteristic value
            omega_hat(i) = abs(poles(ii)); % Absolute value to the characteristic value
            zeta_hat(i) = - real(poles(ii))/omega_hat(i); % -real(lambda)/abs(lambda) (Normalized real part)
            phi(:,i) = modes(:,ii); % Mode shape
        end
        [omega_hat_sort, I] = sort(omega_hat, 'ascend');
        zeta_hat_sort = zeta_hat(I);
        if nargin == 2 % Default output of modes
            modes_s = phi(:,I);                  % Nc*Nc/2 sorted modes (no conjugate pair)
            mp = [omega_hat_sort; zeta_hat_sort]; % modal parameter matrix
            % Return poles in two columns of complex eigenvalues
            poles_s = mp2poles(mp);
        elseif isequal(lower(type), 'jordanpair')
            % Nc*2Nc complex mode sorted according to the oscillation freq
            % First Nc*Nc are the modes related to displacement
            % Second Nc*Nc are the conjugate modes
            modes_s = [phi(1:Nc/2,I), conj(phi(1:Nc/2,I))];
            mp = [omega_hat_sort; zeta_hat_sort];
            poles_s = mp2poles(mp);
            % 2Nc*Nc complex eigenvalues
            J = [diag(poles_s(:,1)), zeros(Nm,Nm);...
                zeros(Nm,Nm), diag(poles_s(:,2))];
            % Reassign the poles with the eigenvalue matrix J
            poles_s = J;
        elseif isequal(lower(type), 'unique')
            % Return only the Nm*Nm complex mode shapes which correspond to the
            % modes of the displacement response
            modes_s = phi(1:Nc,I);
            % Return the unique complex eigenvalues by rejecting the conjugate
            % pairs, Nm*1
            mp = [omega_hat_sort; zeta_hat_sort];
            poles_s = mp2poles(mp);
            poles_s = poles_s(:,1);
        elseif isequal(lower(type), 'uniquedisp')
            % Return only the Nm*Nm complex mode shapes which correspond to the
            % modes of the displacement response
            modes_s = phi(1:Nc/2,I);
            % Return the unique complex eigenvalues by rejecting the conjugate
            % pairs, Nm*1
            mp = [omega_hat_sort; zeta_hat_sort];
            poles_s = mp2poles(mp);
            poles_s = poles_s(:,1);
        end
    end
end

