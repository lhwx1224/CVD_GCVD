function plotDampingErrorStat(mps_stat_mean, markers, colors, alphas, legends, mps_stat_std)
% plotMps visualize the system modal parameters in terms of the natural
% frequencies and damping ratios.
% syntax: plotMps(mps_hat, mps_true, markers, colors, alphas, legends)
% input: mps_hat - cell, in which a set of modal parameter matrices are
%                        stored
%        mps_true - array, in which the true modal parameter matrix is
%                          stored
%        markers - cell, in which markers of each batch of poles are
%        presented
%        colors - cell, in which colors of each batch of poles are used
%        alpha - cell, in which alpha/transparency value of the poles are
%        used
%        legends - cell, a seores of strings that describe the poles
% output: a plot that illustrates the estimated parameter errors.
% Copyright by Hewenxuan Li, 2023 @ Cornell
DefaultFontSize = 10;
set(0, 'defaultAxesFontSize',DefaultFontSize)
% figure
if nargin < 2
    markers = {'o','x','x','x','x'};
    colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
    alphas = [1,1,1,1,1];
    legends = {'Pole1','Pole2','Pole3','Pole4','Pole5'};
end

if nargin > 1 && nargin < 5
    error('All five inputs are needed!')
end

if nargin < 6
    plottype = 'line';
else
    plottype = 'errorbar';
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
    legends = {'Pole1','Pole2','Pole3','Pole4','Pole5'};
end

N = length(mps_stat_mean);

clf
if isequal(lower(plottype), 'line')

for i = 1:N
    mps_hat_i = mps_stat_mean{i}; % Assign the modal parameters
    [m_,n_] = size(mps_hat_i);
    p1(i) = plot(mps_hat_i(3,:),markers{i},'Color',colors{i},'markersize',5,'linewidth',1.5);
    hold on
    p1(i).Color = [p1(i).Color, alphas(i)];
end
else
for i = 1:N
    mps_hat_i = mps_stat_mean{i}; % Assign the modal parameters
    mpsstd_i = mps_stat_std{i};
    [m_,n_] = size(mps_hat_i);
    x = 1:1:n_;
    % Plot the errorbar without markers
    p1(i) = errorbar(x, mps_hat_i(3,:), mpsstd_i(3,:), ...
        ...
        'Color',colors{i}, ...
        'linewidth',1.5);
    hold on
    set([p1(i).Bar, p1(i).Line], 'Colortype','truecoloralpha', 'ColorData', [p1(i).Line.ColorData(1:3); 255*alphas(i)])
    % Impose markers on top of the errorbar using scatter
    s1(i) = scatter(x, mps_hat_i(3,:), 20, hex2rgb(colors(i)), 'Marker',markers{i}, 'LineWidth',1.5, 'MarkerEdgeAlpha',alphas(i));
end    
end
% print('-dpng','-r600',strjoin({FigurePath,['Cmodes_',num2str(n),'dof_mps_err_comp_s-vd_gcvd.png']},filesep))
leg = legend([s1, p1(1)],legends,'Location','Northwest','interpreter','latex','numcolumns',3, 'location', 'Southwest','Fontsize',DefaultFontSize-1);
leg.ItemTokenSize = [13,13]; 
ylabel('$e_\zeta$','Interpreter','latex','Fontsize',DefaultFontSize)
xlabel('Mode Index','Interpreter','latex','Fontsize',DefaultFontSize)
% yticks([1e-4,1e-2,1e0])
% yticklabels({'$10^{-4}$','$10^{-2}$','$10^{0}$'})
if n_ <= 10
    xticks(1:n_+2)
    nums = string(1:n_);
    xticklabels({nums{:}})
else
    xticks(0:dm:n_+2)
    nums = string(1:dm:n_);
    xticklabels({'',nums{:},''})
end
xlim([0.9 10.1])
set(gca,'fontsize',DefaultFontSize)
set(gcf,'Paperposition',[0 0 3 3*0.5])