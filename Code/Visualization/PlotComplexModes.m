function PlotComplexModes(modesn)

% Load environment
ProjName = 'CVD_GCVD';
folders = split(pwd,filesep);
for i = 1:length(folders)
    if isequal(folders{i}, ProjName)
        flag = i;
        break
    end
end
ProjectPath = strjoin(folders(1:flag),filesep);
FigurePath = strjoin({ProjectPath,'Figures','MDOF'}, filesep);
DataPath = strjoin({ProjectPath,'Data','MDOF'}, filesep);
CodePath = strjoin({ProjectPath,'Code'}, filesep);
addpath(CodePath)
addpath(strjoin({ProjectPath, 'OS'},filesep))
addpath(strjoin({CodePath, 'Visualization'},filesep))
addpath(strjoin({CodePath, 'Metrics'},filesep))
checkdir(FigurePath);
[Markers, LineStyles, MyColors] = mystyle();
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultTextInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
% --------------------------------------------
% Begin the 3D complex mode plot
[loc, n] = size(modesn);
loc = 1:loc;
np = 3;
for j = 1:np
    if np > 7
        p{j}=plot3(loc, real(modesn(:,j)), imag(modesn(:,j)),'o-','color',MyColors.defaultcolor{1});
    else
        p{j}=plot3(loc, real(modesn(:,j)), imag(modesn(:,j)),'o-','color',MyColors.defaultcolor{j});
    end
    hold on
    for i = 2:n+1
        if n <= 5
            if np > 7
                plot3([loc(i) loc(i)], [0 real(modesn(i,j))], [0 imag(modesn(i,j))],':','color',MyColors.defaultcolor{1})
            else
                plot3([loc(i) loc(i)], [0 real(modesn(i,j))], [0 imag(modesn(i,j))],':','color',MyColors.defaultcolor{j})
            end
        end
    end
end
plot3(loc,zeros(size(loc)),zeros(size(loc)),'k:')
pbaspect([2 1 1])
grid on
dm = round(n/6);
if n <= 10
    xticks(1:size(modesn,2)+2)
    nums = string(1:n);
    xticklabels({'',nums{:},''})
else
    xticks(0:dm:size(modesn,2)+2)
    nums = string(1:dm:n);
    xticklabels({'',nums{:},''})
end
xlim([1 n+2])
if n > 15
    zlim([-0.35 0.35])
end
xlabel('Spatial Location')
ylabel('$\Re \phi$')
zlabel('$\Im \phi$')
plots = [];
texts = cell(np,1);
for i = 1:np
plots = [plots p{i}];
texts{i} = ['$\phi_{',num2str(i),'n}$'];
end
lgd = legend(plots, texts ,'location','southwest','NumColumns',3);
title('Complex Modes')
set(gcf,'papersize',[6 4])
set(gcf,'paperposition',[0 0 6 4])
% ------------------------------------------------------------
% Plot the modes in the side view (in the complex plane)
figure
for j = 1:np
    if np > 7
        p{j}=plot3(loc, real(modesn(:,j)), imag(modesn(:,j)),'o-','color',MyColors.defaultcolor{1});
    else
        p{j}=plot3(loc, real(modesn(:,j)), imag(modesn(:,j)),'o-','color',MyColors.defaultcolor{j});
    end
    hold on
    for i = 2:n+1
        if np > 7
            plot3([loc(i) loc(i)], [0 real(modesn(i,j))], [0 imag(modesn(i,j))],':','color',MyColors.defaultcolor{1},'LineWidth',0.5)
        else
            plot3([loc(i) loc(i)], [0 real(modesn(i,j))], [0 imag(modesn(i,j))],':','color',MyColors.defaultcolor{j},'LineWidth',0.5)
        end
    end
end
% ylim([rmin+0.1*rmin, rmax+0.1*rmax])
% if n > 15
%     zlim([-0.35 0.35])
% else
%     zlim([imin+0.1*imin, imax+0.1*imax])
% end
ylabel('$\Re \phi$')
zlabel('$\Im \phi$')
title('Complex Plane View')
view([90, 0])
pbaspect([1 1 1])
box on

end