function plotMps(mps_hat, mps_true, markers, colors, alphas, legends)
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

set(0, 'defaultAxesFontSize',12)
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

N = length(mps_hat);

clf
for i = 1:N
    mps_hat_i = mps_hat{i}; % Assign the modal parameters
    % The first plot shows the error in natural freuqency
    subplot(211)
    p1(i) = plot(abs(mps_true(2,1:end)-mps_hat_i(2,:))./mps_true(2,1:end),markers{i},'Color',colors{i},'markersize',6,'linewidth',1.5);
    hold on
    p1(i).Color = [p1(i).Color, alphas(i)];
    % The second plot shows the error in modal damping ratios
    subplot(212)
    p2(i) = plot(abs(mps_true(3,1:end)-mps_hat_i(3,:))./mps_true(3,1:end),markers{i},'Color',colors{i},'markersize',6,'linewidth',1.5);
    hold on
    p2(i).Color = [p2(i).Color, alphas(i)];
end
% print('-dpng','-r600',strjoin({FigurePath,['Cmodes_',num2str(n),'dof_mps_err_comp_s-vd_gcvd.png']},filesep))

subplot(211)
ylabel('$e_f$','Interpreter','latex')
legend(legends,'Location','Northwest','interpreter','latex','numcolumns',2, 'location', 'Southwest')
xlim([0.9 10.1])
subplot(212)
xlabel('Mode Index','Interpreter','latex')
ylabel('$e_\zeta$','Interpreter','latex')
xlim([0.9 10.1])
set(gcf,'paperposition',[0,0,6,4])