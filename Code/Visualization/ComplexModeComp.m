function ComplexModeComp(Phi_hat, Phi, mode_indx, views, dmethod)
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
if nargin < 3
    mode_indx = [1,2,3];
    views = 'default';
    dmethod = '';
elseif nargin < 4
    views = 'default';
    dmethod = '';
elseif nargin < 5
    dmethod = '';
end
[Markers, LineStyles, MyColors] = mystyle();
n = size(Phi, 2);
% match the phase to the actual modes
Phi_hat_pm = ModePhaseMatch(Phi_hat, Phi, n);

% Visualization
loc = 1:n+2;

figure,clf
count = 1;
for j = 1:length(mode_indx)
    p1(count)=plot3(loc, real(Phi(:,mode_indx(j))), imag(Phi(:,mode_indx(j))),'o-','color',MyColors.defaultcolor{j});
    p1(count).Color = [p1(count).Color, 0.5];
    p1(count).LineWidth = 2;
    hold on
    p2(count)=plot3(loc, Phi_hat_pm{mode_indx(j)}(:,1), Phi_hat_pm{mode_indx(j)}(:,2),'*--','color',MyColors.defaultcolor{j});
    count = count + 1;
% The vector line plot that connects the zero position and the nodal value
%     for i = 2:n
%         plot3([loc(i) loc(i)], [0 real(Phi(i,j))], [0 imag(Phi(i,j))],':','color',MyColors.defaultcolor{j})
%         hold on
%         plot3([loc(i) loc(i)], [0 Phi_hat_pm{j}(i,1)], [0 Phi_hat_pm{j}(i,2)],'--','color',MyColors.defaultcolor{j+1})
%     end
end
plot3(loc,zeros(size(loc)),zeros(size(loc)),'k:')
pbaspect([2.5 1 1])
grid on
xlabel('Spatial Location', 'interpreter','latex')
ylabel('$\Re \phi$', 'interpreter', 'latex')
zlabel('$\Im \phi$', 'interpreter', 'latex')
plots = [];
texts = cell(length(mode_indx)*2,1);
for i = 1:length(mode_indx)
plots = [plots p1(i)];
plots = [plots p2(i)];
texts{2*i-1} = ['$\phi_{',num2str(mode_indx(i)),'n}$'];
texts{2*i} = ['$\hat\phi_{',num2str(mode_indx(i)),'n\mathrm{',dmethod,'}}$'];
end
lgd = legend(plots, texts ,'location','North','NumColumns',length(mode_indx), 'interpreter', 'latex');

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