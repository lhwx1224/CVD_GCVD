function ComplexModePlot_yy(Phi, mode_indx)
% ComplexModePlot_yy - plots the complex mode shapes in terms of their real
% and imagniary parts.
% syntax: ComplexModePlot_yy(Phi, mode_indx)
% input: Phi - array, complex mode shapes matrix
%        mode_indx - the index of the modes to be plotted
% output: Plot of modes overlay on a yyaxis plot
%
% Hewenxuan Li @ Cornell 04-15-2023
if nargin < 2
    mode_indx = [1,2,3];
end

if isreal(Phi)
    error('Modes have to be complex-valued')
end
n = size(Phi,2);
np = length(mode_indx); % Number of plotted modes
if np > n
    np = n;
end
% Get the min and max range of the modes
rmax = max(real(Phi(:,mode_indx)),[],'all');
rmin = min(real(Phi(:,mode_indx)),[],'all');
imax = max(imag(Phi(:,mode_indx)),[],'all');
imin = min(imag(Phi(:,mode_indx)),[],'all');

rran = max(rmax, abs(rmin));
iran = max(imax, abs(imin));

plots = zeros(np*2,1);
legends = cell(np*2,1);
figure, clf
yyaxis left
for i = 1:np
    pl(i) = plot(real(Phi(:,mode_indx(i))));
    hold on
    plots(2*i-1) = pl(i);
    legends{2*i-1} = ['$\Re \phi_', num2str(i), '$'];
end
xticks([1, 2, 3, 4])
xticklabels({'','','',''})
ylabel('$\Re \phi$')
ylim(1.2*[-rran rran])

yyaxis right
for i = 1:np
    pr(i) = plot(imag(Phi(:,mode_indx(i))));
    plots(2*i) = pr(i);
    legends{2*i} = ['$\Im \phi_', num2str(i), '$'];
end
ylabel('$\Im \phi$')
pbaspect([2 1 1])
xticks(1:n+2)
nums = string(1:n);
xticklabels({'',nums{:},''})
xlim([1 n+2])
ylim(1.2*[-iran iran])
xlabel('Spatial Location')
legend(plots, legends,'numcolumns',np)
set(gcf,'papersize',[4 2])
set(gcf,'paperposition',[0 0 4 2])
set(gcf,'renderer','painters')