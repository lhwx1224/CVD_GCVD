function ComplexModeComp_(Phi_hat_list, Phi, mode_indx, views, dmethod, markers, colors)
% ComplexModeComp
% syntax: Phi_hat - complex estimated modes
%         Phi - complex true modes
%         mode_indx - the index of the modes to be plotted
%         views - string, selected views for visualization
%                 'default' - MATLAB default 3D view angle
%                 'real' - look through the real part of the modes
%                 'imag' - look through the imaginary part of the modes
%                 'complexplane' - look through the complex plane direction
% Hewenxuan Li @ Cornell 04-10-2023
DefaultFontSize = 10;
set(0, 'defaultAxesFontSize',DefaultFontSize)
if nargin < 3
    mode_indx = [1,2,3];
    views = 'default';
    dmethod = '';
    markers = {'^','v','d','s','x'};
elseif nargin < 4
    views = 'default';
    dmethod = '';
    markers = {'^','v','d','s','x'};
elseif nargin < 5
    markers = {'^','v','d','s','x'};
    dmethod = '';
elseif nargin < 6
    markers = {'^','v','d','s','x'};
end
[Markers, LineStyles, MyColors] = mystyle();
% Placeholder for the plots
plots = [];
N = length(Phi_hat_list);
n = size(Phi, 2);
Ni = length(mode_indx);
% Create a new figure
figure(1),clf
loc = 1:n+2;
% The very first plot is the ground truth
count = 1;
for j = 1:Ni
    if length(colors) == 0
        p{1}(count)=plot3(loc, real(Phi(:,mode_indx(j))), imag(Phi(:,mode_indx(j))),'o-','color',MyColors.defaultcolor{j});
    else
        p{1}(count)=plot3(loc, real(Phi(:,mode_indx(j))), imag(Phi(:,mode_indx(j))),'o-','color',colors{j});
    end
    p{1}(count).Color = [p{1}(count).Color, 0.5];
    p{1}(count).LineWidth = 2;
    count = count + 1;
    hold on
end
for k = 1:N
% match the phase to the actual modes
Phi_hat = Phi_hat_list{k};
Phi_hat_pm = ModePhaseMatch(Phi_hat, Phi, n);
% Plot the estimated modes
count = 1;
for j = 1:Ni
    if length(colors) == 0
    p{k+1}(count)=plot3(loc, Phi_hat_pm{mode_indx(j)}(:,1), Phi_hat_pm{mode_indx(j)}(:,2),[markers{k},'--'],'color',MyColors.defaultcolor{j});
    else
    p{k+1}(count)=plot3(loc, Phi_hat_pm{mode_indx(j)}(:,1), Phi_hat_pm{mode_indx(j)}(:,2),[markers{k},'--'],'color',colors{j});
    end
    p{k+1}(count).Color = [p{k+1}(count).Color 0.7];
    count = count + 1;
end
end
texts = cell(Ni*(N+1),1);
for i = 1:Ni
    for k = 1:N+1
        plots = [plots p{k}(i)];
        if k == 1
            texts{(N+1)*(i-1)+k} = ['$\phi_{',num2str(mode_indx(i)),'n}$'];
        else
            texts{(N+1)*(i-1)+k} = ['$\hat\phi_{',num2str(mode_indx(i)),'n\mathrm{,\,',dmethod{k-1},'}}$'];
        end
    end
end
% Plot the equilibrium line
plot3(loc,zeros(size(loc)),zeros(size(loc)),'k:')
lgd = legend(plots, texts ,'location','North','NumColumns',Ni, 'interpreter', 'latex','FontSize',8);
leg.ItemTokenSize = [13,13]; 
pbaspect([2.5 1 1])
grid on
xlabel('Spatial Location', 'interpreter','latex')
ylabel('$\Re \phi$', 'interpreter', 'latex')
zlabel('$\Im \phi$', 'interpreter', 'latex')
xticks([2:11])
xticklabels(string(1:10))
xlim([1,12])


% Assign the specified view
if isequal(lower(views), 'default')
    view([-37.5, 30])
elseif isequal(lower(views), 'real')
    view([0, 90])
elseif isequal(lower(views), 'imag')
    view([0, 0])
elseif isequal(lower(views), 'complexplane')
    view([90, 0])
end
set(gcf,'papersize',[5 3])
set(gcf,'paperposition',[0 0 5 3])
set(gcf,'renderer','painters')