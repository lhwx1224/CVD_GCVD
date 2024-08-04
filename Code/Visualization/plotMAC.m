function plotMAC(MACs, markers, colors, alphas, legends, MACs_std)
% plotMAC visualize the calculated modal assurances for each estimated set
% of modes. The modal assurances are obtained from the modal assurance
% criterion (MAC), please see MAC.mat for more details.
% syntax: plotMAC(MACs, markers, colors, alphas, legends)
% input: MACs - cell, in which a set of MAC matrices are stored
%        markers - cell, in which markers of each batch of poles are
%        presented
%        colors - cell, in which colors of each batch of poles are used
%        alpha - cell, in which alpha/transparency value of the poles are
%        used
%        legends - cell, a seores of strings that describe the poles
% output: a plot that illustrates the MAC values.
% Copyright by Hewenxuan Li, 2023 @ Cornell

DefaultFontSize = 10;
set(0, 'defaultAxesFontSize',DefaultFontSize)
% figure
if nargin < 3
    markers = {'o','x','x','x','x'};
    colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
    alphas = [1,1,1,1,1];
    legends = {'Method1','Method2','Method3','Method4','Method5'};
end

if nargin < 6
    plottype = 'line';
else
    plottype = 'errorbar';
end

if nargin > 2 && nargin < 5
    error('All five inputs are needed!')
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


N = length(MACs);

% clf
% for i = 1:N
%     mac_i = diag(MACs{i}); % Assign the modal assurance criteria 
%     p1(i) = plot(mac_i,markers{i},'Color',colors{i},'markersize',5,'linewidth',1.5);
%     hold on
%     p1(i).Color = [p1(i).Color, alphas(i)];
% end
if isequal(lower(plottype), 'line')
for i = 1:N
    mac_i = diag(MACs{i}); % Assign the modal parameters
    p1(i) = plot(mac_i,[markers{i},'-'],'Color',colors{i},'markersize',5,'linewidth',1.5);
    hold on
    p1(i).Color = [p1(i).Color, alphas(i)];
end
leg = legend([p1],legends,'Location','Northwest','interpreter','latex','numcolumns',3, 'location', 'SouthEast','Fontsize',DefaultFontSize-1);
else
for i = 1:N
    mac_i = diag(MACs{i}); % Assign the modal parameters
    mpsstd_i = diag(MACs_std{i});
    [m_,n_] = size(mac_i);
    x = 1:1:m_;
    p1(i) = errorbar(x, mac_i, mpsstd_i, ...
        ...
        'Color',colors{i}, ...
        'linewidth',1.5);
    hold on
    set([p1(i).Bar, p1(i).Line], 'Colortype','truecoloralpha', 'ColorData', [p1(i).Line.ColorData(1:3); 255*alphas(i)])
    s1(i) = scatter(x, mac_i, 20, hex2rgb(colors(i)), 'Marker',markers{i}, 'LineWidth',1.5, 'MarkerEdgeAlpha',alphas(i));
end    
leg = legend([s1],legends,'Location','Northwest','interpreter','latex','numcolumns',3, 'location', 'SouthEast','Fontsize',DefaultFontSize-1);
end
leg.ItemTokenSize = [13,13];
xlim([0.9 10.1])
xlabel('Mode Index','Interpreter','latex','FontSize',DefaultFontSize)
ylabel('MAC','Interpreter','latex','FontSize',DefaultFontSize)
set(gca,'Fontsize',DefaultFontSize)
xticks(1:10)
xticklabels(string(1:10))
xlim([0.9 10.1])
% ylim([0 1])
set(gcf,'paperposition',[0,0,4,2])
set(gcf,'Paperposition',[0 0 3 3*0.5])