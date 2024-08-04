function ComplexModePlot(Phi, mode_indx)
% ComplexModePlot - plots the complex mode shapes in terms of their real
% and imagniary parts.
% syntax: ComplexModePlot(Phi, mode_indx)
% input: Phi - array, complex mode shapes matrix
%        mode_indx - the index of the modes to be plotted
% output: Plot of modes overlay
%
% Hewenxuan Li @ Cornell 04-15-2023

[Markers, LineStyles, MyColors] = mystyle();
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultTextInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
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
ran = max(rran, iran);
plots = zeros(np*2,1);
legends = cell(np*2,1);
figure(1),clf
for i = 1:np
    if np > 7
        pr(i)=plot(real(Phi(:,i)),'color',MyColors.defaultcolor{1},'LineStyle',LineStyles{1});
        hold on
        pi(i)=plot(imag(Phi(:,i)),'color',MyColors.defaultcolor{2},'LineStyle',LineStyles{2});
    else
        pr(i)=plot(real(Phi(:,i)),'color',MyColors.defaultcolor{i},'LineStyle',LineStyles{1});
        hold on
        pi(i)=plot(imag(Phi(:,i)),'color',MyColors.defaultcolor{i},'LineStyle',LineStyles{2});
    end
    plots(2*i-1) = pr(i);
    plots(2*i) = pi(i);
    legends{2*i-1} = ['$\Re \phi_', num2str(i), '$'];
    legends{2*i} = ['$\Im \phi_', num2str(i), '$'];
end
xticks(1:n+2)
nums = string(1:n);
xticklabels({'',nums{:},''})
xlabel('Spatial Location')
ylabel('Magnitude')
pbaspect([2 1 1])
legend(plots, legends, 'numcolumns', np)
xlim([1 n+2])
ylim(1.5*[-ran ran])
set(gcf,'papersize',[4 2])
set(gcf,'paperposition',[0 0 4 2])
set(gcf,'renderer','painters')