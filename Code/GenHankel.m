function H = GenHankel(Y, nd)
% GenHankel(Y, nr) takes in a data matrix y and convert it into a Hankel
% matrix based on how many rows are desired, nr. The general Henkel matrix
% can be constructed by imposing more constaints, which are (1) the delay
% time, tau, used in the Hankel matrix construction; this number is defalt
% set to tau = 1; and (2) the resampling distance used, r. 
%
% Syntax: H = GenHankel(y, nr)
% 
% Input: y - 2d array, whose rows are the channel indices used in the data
%            acquisition.
%        nd - integer, number of delays used for constructing the Hankel
%             structure.
%
% Copyright by Hewenxuan Li, April 2021, for Modal Identification Project
%
% Created on April 20, 2021
% Last modified:

if size(Y, 1) > size(Y,2)
    Y = Y';
end

[n, m] = size(Y);

% The Hankel matrix is composed of nd+1 blocks with m-nd time columns
H = zeros((nd+1)*n, m-nd);

% Fill in the Hankel matrix
for i=1:nd + 1
    H(n*(i-1)+1:i*n,:) = Y(:,i:m - (nd - i + 1)); % Assign the delayed vectors to Hankel Structure
end