function mp = poles2mp(poles, varargin)
%POLES2MP Converting poles to modal parameters
%   MP = POLES2MP(S) returns the modal parameter matrix, mp, whose first
%   row contains the angular natural frequency of the identified modes, the
%   second row contains the identified natural frequencies in Hz, and the
%   third row contains the damping ratios, zeta.
%
%   Syntax: MP = POLES2MP(S)
%   
%   Input: S - complex velued vector, estimated/analytical poles 
%
%   Output: MP - real valued matrix, 3*NM, estimated/analytical modal
%   parameters: 
%       row1 contains the angular natural frequencies, 
%       row2 contains the natural frequency in Hz, and 
%       row3 contains the damping ratios, zeta.
%
% Copyright by Hewenxuan Li, April 2021
% References: [1] Fundamentals of Vibrations, 2000, Meirovitch (see p.88 Eq. 2.25)

% Check if the input is a matrix or a vector, if a matrix, turn into vector
if size(poles,1) ~=1 || size(poles,2) ~= 1
    poles = diag(poles);
end
N = length(poles);
mp = zeros(3, N);                         % Modal parameter estimate

if ~isreal(poles)
    for i=1:N
        omegan = abs(poles(i));              % angular natural frequency
        mp(1,i) = omegan;
        mp(2,i) = omegan/2/pi;           % frequency in Hz
        mp(3,i) = -real(poles(i))/omegan;    % zeta
    end
    
    if nargin > 1
        if isequal(lower(varargin{1}), 'unique')
            mp = mp(:,1:2:end);
        end
    end
else % If the poles are purely real - undamped system
    for i=1:N
        omegan = sqrt(poles(i));         % angular natural frequency
        mp(1,i) = omegan;
        mp(2,i) = omegan/2/pi;           % frequency in Hz
        mp(3,i) = 0;    % zeta
    end
end