function [Phi_hat_pm, Is, I, Sij, Pij, Sscore, Pscore] = genModePhaseMatch(Phi_hat, Phi, n, method)
% genModePhaseMatch matches the phases between the estimated modes
% 'Phi_hat' and the true modes 'Phi' for 'n' number of modes; the
% phase-matched mode shapes are denoted as 'Phi_hat_pm'. The method is
% based on a psudo-inverse as a rotational transformation to the estimated
% mode shape. 
% syntax:  [Phi_hat_pm, Is, I, Sij, Pij, Sscore, Pscore] = genModePhaseMatch(Phi_hat, Phi, n)
% input: Phi_hat - Estimated mode shapes
%        Phi - true mode shapes
%        n - number of modes to match
% ouput: Phi_hat_pm - phase-matched estimated mode shapes, sorted according
%        to the ascending true mode index, using the Is index.
%        Is - sorting index that used to sort the esitmated mode shapes
%        I - index that indicates with respect to which mode the estimated
%        modes are the most similar
%        Sij - index matrix whose element s_{i,j} is the index of the
%        estimated modes that are most relavent to the ith true mode, the
%        index j indicates the strength of similarity in descending orders
%        Pij - index matrix whose element p_{i,j} is the index of the true
%        modes that are most relavant to the ith esitimated modes, the
%        index j indicates the strength of similarity in descending orders.
%        Sscore - similarity score matrix whose i,jth element represents
%        the similarity between the ith true mode and the jth estimated
%        mode.
%        Pscore - similarity score matrix whose i,jth element represents
%        the similarity between the ith estimated mode and the jth true
%        mode.
% Hewenxuan Li 10-08-2023 @ Cornell

% Define the phase-matched mode vectors
if nargin < 3 
    n = size(Phi,1);
    method = 'ReIm';
elseif nargin < 4
    method = 'ReIm';
end

[m_hat, K] = size(Phi_hat);
[~, L] = size(Phi);

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
Phi_hat_pm_cell = cell(K, L);
Pscore = zeros(K,L);
% Evaulate the phase match and obtain the matching matrix Pscore (phase-match score)
for k = 1:K
    % ------------------------ Force modes to be complex -------------
    if isreal(Phi_hat(:,k))
        Phi_hat(:,k) = Phi_hat(:,k) + 1i*1e-6*randn(size(Phi_hat(:,k)));
    end
    % ---------------------------------------------------------------
    for l = 1:L
    PT = ([real(Phi_hat(:,k)), imag(Phi_hat(:,k))]'*[real(Phi_hat(:,k)), imag(Phi_hat(:,k))])\...
         [real(Phi_hat(:,k)), imag(Phi_hat(:,k))]'*[real(Phi(:,l)) imag(Phi(:,l))]; % Phase matching transform, psudo-inverse based
    rim = [real(Phi_hat(:,k)) imag(Phi_hat(:,k))]*PT;                              % 
    Phi_hat_pm_cell{k,l} = rim(:,1) + 1i*rim(:,2);
    Pscore(k,l) = Phi_hat_pm_cell{k,l}'*Phi(:,l)/(norm(Phi_hat_pm_cell{k,l})*norm(Phi(:,l)));
%     if sum(isnan(Phi_hat_pm{k}),'all') ~= 0
%         Phi_hat_pm{k} = Phi_hat{k};
%     end
    end
end
Pscore = real(Pscore); % Convert the matrix to real-valued
Sscore = Pscore';      % Sorting score matrix
% Find the matching index and sort the Phi_hat
I = zeros(K,1);            % initialize the best-matched index of true mode to the estimated mode
Is = zeros(L,1);           % initialize the bset-sorting index from the estimated mode to the true mode
Pij = zeros(size(Pscore)); % Initialize the sorted index simlarity matrix
Sij = zeros(size(Sscore));
Phi_hat_pm = zeros(m_hat, min(K,L)); % Initialize the sorted estimated modes, size depends on the mininum value between the two matrices

for k = 1:K
    [~, indx] = sort(Pscore(k,:),"descend"); % Sort the Pscore in descending order
    Pij(k,:) = indx;  % Assign the sorted index (the similarity between the phi_hat_i and phi_j)
    I(k) = indx(1);   % Assign the best-matched index to phi_hat_i
end
for l = 1:L
    [~, indxs] = sort(Sscore(l,:), "descend"); % Sort the Sscore in descending order
    Sij(l,:) = indxs;
    Is(l) = indxs(1);
end
% Construct the sorted modes
for k = 1:min(K,L)
    % Find the best-phase-matched candidate according to Is and assign to the new mode shape matrix
    Phi_hat_pm(:,k) = Phi_hat_pm_cell{Is(k),k}; 
end
end
