function [C, ModeIDIndex, SortIndex, varargout] = MAC(Phi_hat, Phi_a)
% MAC returns the modal assurance criterion matrix C and the sorted index
% for the analytical mode that is corresponding to the estimated mode
% shape. If no output variables specified, the matrix plot and a 3-d bar
% plot will be returned, indicating the quality of the modal
% identification. 
% Note that for complex modes, the genModePhaseMatch will first find the
% best-phase-matched modes and use the best-matched modes to do MAC with
% repsect to the true modes. Hence, an equivalent raw MAC value may be
% obtained by calling varargout{1} for Sij as the traditional MAC value.
% 
% Syntax: [C, ModeIDIndex, SortIndex] = MAC(Phi_hat, Phi_a)
%
% Input:    Phi_hat - estimated mode shapes
%             Phi_a - analytical/FEM mode shapes
%
% Output:         C - MAC matrix, correlation between estimated modes and
%                     the reference modes
%       ModeIDIndex - index indicating the closest reference mode shape to
%                     the currently considered esimated modes. If the input
%                     matrix is real, the MAC value will be assiend by
%                     finding the maximum value in the MAC matrix in the
%                     column direction; if the input matrix is complex, the
%                     index is obtained according to I based on phase
%                     matching similarity Pij.
%         SortIndex - index indicating how the estimated modes should be
%                     sorted to follow the reference mode order. If the
%                     input matrix is real, the MAC value will be assiend
%                     by finding the maximum value in the MAC matrix in the
%                     row direction; if the input matrix is complex, the
%                     index is obtained according to Is based on phase
%                     sorting similarity Sij.
%
% Copyright by Hewenxuan Li, 2021 
% Last Modified: 
% v1.0 by Hewenxuan Li, May 2023 @ Cornell - extended to complex modes
% v2.0 by Hewenxuan Li, Oct 2023 @ Cornell - incorporated generalized phase
%                                            matching algorithm
% References: [1] Modal Assurance Criterion - Pastor et. al., 2012

[~, n] = size(Phi_hat);

if ~isreal(Phi_hat)
%     Phi_hat_ri = ReImModes(Phi_hat);
%     Phi_a_ri = ReImModes(Phi_a);
%     Phi_hat_mp_ri = ModePhaseMatch(Phi_hat_ri, Phi_a, size(Phi_hat_ri,1),'ReIm');
    Phi_hat_temp = Phi_hat;
    [Phi_hat, Is, I, Sij, Pij] = genModePhaseMatch(Phi_hat, Phi_a, size(Phi_hat,2),'complex');
    % Check if all modes are properly matched
    indx = find(sum(isnan(Phi_hat),1),n+2);
    if ~isempty(indx)
        warning(['Phase-mathed modes failed! mode index:', num2str(indx), ' are detected to be NaNs!'])
        warning('NaN valued modes are replaced by the original mode estimates!')
        Phi_hat(:,indx) = Phi_hat_temp(:,indx);
    end
end

% Calculate the vector norm to each mode for both analytical and estimated.
norm_hat = vecnorm(Phi_hat);
norm_a = vecnorm(Phi_a);

% Make them squared
norm_hat = norm_hat.^2; % First term of the denominator
norm_a = norm_a.^2;     % Second term of the denominator

% Span the denominator matrix
D = norm_hat'*norm_a;

% Calculate the cross-covariance of the estimated modes and the analytical
% ones
N = Phi_hat'*Phi_a;

% Make it squared
if isreal(N)
    N = N.^2;
else
    N = abs(N).^2;
end

% MAC Matrix
C = N./D;

if isreal(Phi_hat)
% Index indicating the closest reference mode to the estimated modes
[~, ModeIDIndex] = max(C, [], 2);
% Index indicating the closest estimated mode to the reference modes
[~, SortIndex] = max(C, [], 1);
varargout = {};
else
    ModeIDIndex = I;
    SortIndex = Is;
    varargout = {Sij, Pij};
end

if nargout == 0
    figure
    subplot(121)
    imagesc(C)
    colorbar
    % DEFINE THE VIRIDIS COLOR MAP
    color1 = '#FDE725';
    color2 = '#7AD151';
    color3 = '#22A884';
    color4 = '#2A788E';
    color5 = '#414487';
    color6 = '#440154';  
    map = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6});
    colormap(map)
    pbaspect([1 1 1])
    xlabel('Index of $\hat{\Phi}$','interpreter','latex')
    ylabel('Index of ${\Phi}$','interpreter','latex')
    subplot(122)
    h = bar3(C);
    % -------------- CODE FOR INTERPOLATED BAR COLOR ---------------------
    shading interp
    for i = 1:length(h)
        zdata = get(h(i),'Zdata');
        set(h(i),'Cdata',zdata)
        set(h,'EdgeColor','k')
    end
    set(h,'EdgeColor', 'none','LineWidth', 0.1);
% -------------------- CODE FOR SOLID BAR COLOR --------------------------
%     numBars = size(C,1);
%     numSets = size(C,2);
%     for i = 1:numSets
%         zdata = ones(6*numBars,4);
%         k = 1;
%         for j = 0:6:(6*numBars-6)
%             zdata(j+1:j+6,:) = C(k,i);
%             k = k+1;
%         end
%         set(h(i),'Cdata',zdata)
%     end
%     
% end
xlabel('Index of $\hat{\Phi}$','interpreter','latex')
ylabel('Index of ${\Phi}$','interpreter','latex')
zlabel('MAC')
pbaspect([1 1 1])
axis tight
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])
set(gcf, 'renderer', 'opengl')
end