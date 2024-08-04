function ReImPhi = ReImModes(Phi, SystemCategory)
% ReImModes takes the provided modes Phi and decompose the modes into their
% real and imaginary parts Re(Phi) and Im(Phi), and outputs a cell array
% that contains the real and imagniary parts as {Re(Phi), Im(Phi)}.
%
% syntax: ReImPhi = ReImModes(Phi, SystemCategory)
% input: Phi - comlex-valued array, complex mode shapes
%        SystemCategory - string, category of the system
%                         if not provided -> 'generic' will be assigned
% output: ReImPhi - cell, the real-imaginary representation of the complex
%                   mode.
%
% Created by Hewenxuan Li @ Cornell, 04-15-2023

% Check the system category
if nargin<2
    SystemCategory = 'generic';
end

% Get the number of modes
n = size(Phi, 2);

% Append boundary conditions if SystemCategory is not generic
if isequal(lower(SystemCategory), 'mdof_ff')
    Phi = [zeros(1,n); Phi; zeros(1,n)];
end

% Initialize the real-imaginary pair as cell an array
ReImPhi = cell(n,1);

% Assign the real and imag part to the cell array in a loop
for i = 1:n
    ReImPhi{i} = [real(Phi(:,i)), imag(Phi(:,i))];
end
