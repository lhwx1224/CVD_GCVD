function PlotModeError(ephi, metric, markers, colors, alphas, legends)
% PlotPoles visualize the modal identification error 
% syntax: PlotModeError(poles, markers, colors, alphas, legends)
% input: ephi - cell, in which a vector of complex poles are stored
%        metric - string, the method used to quantify the error
%        markers - cell, in which markers of each batch of poles are
%        presented
%        colors - cell, in which colors of each batch of poles are used
%        alpha - cell, in which alpha/transparency value of the poles are
%        used
%        legends - cell, a seores of strings that describe the poles
% output: a plot that shows the input poles in the complex plane
% Copyright by Hewenxuan Li, 2023 @ Cornell

DefaultFontSize = 9;
set(0, 'defaultAxesFontSize',DefaultFontSize)
% figure
if nargin < 3
    markers = {'o','x','x','x','x'};
    colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
    alphas = [1,1,1,1,1];
    legends = {'Method1','Method2','Method3','Method4','Method5'};
end

if nargin > 2 && nargin < 6
    error('All six inputs are needed!')
end

if sum(size(markers)) == 0
    markers = {'o','x','x','x','x'};
end

if sum(size(colors)) == 0
    colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
end

if sum(size(alphas)) == 0
    alphas = [1,1,1,1,1];
end

if sum(size(legends)) == 0
    legends = {'Method1','Method2','Method3','Method4','Method5'};
end


[m, n] = size(ephi);
N = max(m,n);
for i = 1:N
    ephi_temp = ephi{i};
    [mi, ni] = size(ephi_temp);
    n = max(mi,ni);
    if ni ~= 1
        ephi_temp = diag(ephi_temp);
    end
        p(i) = plot(1:n, ephi_temp,['-',markers{i}],'Color',colors{i});
        p(i).Color = [p(i).Color alphas(i)];
    if i == 1
        hold on
    end
end
GCA = gca;
plot([0, 0],[GCA.YLim(1), GCA.YLim(2)],'Linewidth',1,'Color','k','LineStyle',':')
plot([GCA.XLim(1), GCA.XLim(2)],[0,0],'Linewidth',1,'Color','k','LineStyle',':')
set(gcf, 'renderer', 'painters')
pbaspect([2 1 1])
xlabel('Mode Index','interpreter','latex','Fontsize',DefaultFontSize)
if isequal(metric, 'rell2')
    ylabel('Relative $L_2$ Error','interpreter','latex','Fontsize',DefaultFontSize)
elseif isequal(metric, 'rmse')
    ylabel('$RMSE$','Fontsize',DefaultFontSize)
elseif isequal(metric, 'angle')
    ylabel('Subspace Angle ($rad$)','interpreter','latex','Fontsize',DefaultFontSize)
end
set(gca,'FontSize',12)
dm = round(n/6);
if n <= 10
    xticks(1:n+2)
    nums = string(1:n);
    xticklabels({nums{:}})
else
    xticks(0:dm:n+2)
    nums = string(1:dm:n);
    xticklabels({'',nums{:},''})
end
xlim([0.9 n+0.1])
legend(p,legends,'Interpreter','latex','Location','Southeast','NumColumns',2,'Fontsize',DefaultFontSize)
set(gca,'fontsize',DefaultFontSize)
set(gcf,'Paperposition',[0 0 4 2])
box on
end