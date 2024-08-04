function [Phi, Lambda, b] = DMD(Y1,Y2,r)
% DMD executes the dynamic mode decomposition to the given data sets Y1 and
% Y2
[U,Sigma,V] = svd(Y1,'econ');      % Step 1
Ur = U(:,1:r);
Sigmar = Sigma(1:r,1:r);
Vr = V(:,1:r);

Atilde = Ur'*Y2*Vr/Sigmar;    % POD Projected Dynamic Matrix
[W,Lambda] = eig(Atilde);     % Step 3

Phi = Y2*(Vr/Sigmar)*W;       % Step 4
alpha1 = Sigmar*Vr(1,:)';
b = (W*Lambda)\alpha1;