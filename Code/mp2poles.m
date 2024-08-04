function pole = mp2poles(modal_param)
% pole = mp2poles(modal_param) converts a set of modal parameters into
% poles of the system.
%
% syntax: pole = mp2poles(modal_param)
%
% input: modal_param - 2*Nm array where Nm is the number of modes, the
%                      first row stores the natural frequencies of each
%                      mode and the second row stores the modal damping
%                      ratios
%                                    [ angular natural freq. ]
%                                    [ modal damping ratios  ]
%
% output: pole - Nm*2 array whose elements are the poles of the system; the
%                first column being the poles in the positive side of the
%                frequency domain, the second column stores the conjugate
%                poles which live in the negative side of the frequency
%                domain
%
%
% Hewenxuan Li, Mar 2021
% Modified on July 2021

if size(modal_param, 1)~=2 && size(modal_param, 2)~=2
    error('Check the size of the input!')
end
Nm = size(modal_param,2);
omega_n = modal_param(1,:);      % Natural frequencies in radian
zeta = modal_param(2,:);         % Modal dampling coefficients
pole = zeros(Nm,2);
for i = 1:Nm
    pole(i,1) = omega_n(i)*(-zeta(i) + 1i*sqrt(1 - zeta(i)^2)); % calculate poles
    pole(i,2) = conj(pole(i,1)); % calcutate the conjugate pair pole
end
end