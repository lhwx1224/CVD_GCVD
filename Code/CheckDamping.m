function dampingtype = CheckDamping(M, C, K, n)
% syntax: dampingtype = CheckDamping(M, C, K, n)
% input: M - mass matrix
%        C - damping matrix
%        K - stiffness matrix
%        n - number of degree-of-freedoms (DOF)
% output: dampingtype - string, type of damping
%                       'p' - proportional damping
%                       'np' - non-proportional damping
% Check if the # of DOF is provided
if nargin < 4
    n = size(M,1);
end
% Check if the system is proportionally damped
if sum(K*inv(M)*C == C*inv(M)*K,"all") == n^2
    disp('System is proportionally damped!')
    dampingtype = 'p';
else
    disp('System is non-proportionally damped!')
    dampingtype = 'np';
end