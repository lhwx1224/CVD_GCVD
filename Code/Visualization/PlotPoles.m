function PlotPoles(poles, markers, colors, alphas, legends)
% PlotPoles visualize the system poles in the complex plane
% syntax: PlotPoles(poles, markers, colors, alphas, legends)
% input: poles - cell, in which a vector of complex poles are stored
%        markers - cell, in which markers of each batch of poles are
%        presented
%        colors - cell, in which colors of each batch of poles are used
%        alpha - cell, in which alpha/transparency value of the poles are
%        used
%        legends - cell, a seores of strings that describe the poles
% output: a plot that shows the input poles in the complex plane
% Copyright by Hewenxuan Li, 2023 @ Cornell
enumerate = 1; % Debugging parameter to show the index of identified poles
set(0, 'defaultAxesFontSize',12)
figure
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


[m, n] = size(poles);
N = max(m,n);
for i = 1:N
    poles_temp = poles{i};
    [mi, ni] = size(poles_temp);
    if ni ~= 1
        poles_temp = diag(poles_temp);
    end
    %     p(i) = plot(poles_temp,'LineStyle','none','Marker',markers{i},'MarkerEdgeColor',colors{i},'MarkerSize',8,'LineWidth',2);
    %     p(i).Color = [p(i).Color, alphas(i)];
    p(i) = scatter(real(poles_temp), imag(poles_temp),40, markers{i},'MarkerEdgeColor',colors{i},'LineWidth',1.5,'MarkerEdgeAlpha',alphas(i));
    if i == 1
        hold on
    end
end

if enumerate == 1
for i = 1:N
    poles_temp = poles{i};
    if i == 1
    ran_val = max(real(poles_temp))-min(real(poles_temp));
    else
        ran_val = [ran_val, max(real(poles_temp))-min(real(poles_temp))];
        if ran_val(i) < ran_val(i-1)
            ran_val(i) = ran_val(i-1);
        end
    end
    [mi, ni] = size(poles_temp);
    if ni ~= 1
        poles_temp = diag(poles_temp);
    end
    [mp, np] = size(poles_temp);
    % for j = 1:mp
    %     text(real(poles_temp(j))+0.1*(5-i)*ran_val(i), imag(poles_temp(j)),num2str(j),'Color',colors{i})
    % end
end

end
GCA = gca;
plot([0, 0],[GCA.YLim(1), GCA.YLim(2)],'Linewidth',1,'Color','k','LineStyle',':')
plot([GCA.XLim(1), GCA.XLim(2)],[0,0],'Linewidth',1,'Color','k','LineStyle',':')
set(gcf, 'renderer', 'painters')
pbaspect([1 1 1])
xlabel('$\Re(\lambda)$','Interpreter','latex','FontSize',14)
ylabel('$\Im(\lambda)$','Interpreter','latex','FontSize',14)
legend(p,legends,'Interpreter','latex','Location','West')
set(gcf,'Paperposition',[0 0 4 4])
box on
end