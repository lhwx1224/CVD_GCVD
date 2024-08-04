function [Markers, LineStyles, Mycolors] = mystyle(options)
% mystyle set up the displaying and plotting environment in MATLAB.
% Syntax: [Markers, LineStyles, Mycolors] = mystyle();
%
% Output: Markers = {'O','+','x','*','s','d','v','^','>','<','p','h'};
%         LineStyles = {'-','--','-.',':'};
%         MyColors.defaultcolor = MATLAB Default
%         MyColors.mycolor -> brighter personal preferences
%         MyColors.mycolor -> duller personal preference
%         MyColors.greys -> greyscale from black to white
%
% Activated setup:
%         1. Figure color set to white
%         2. Default LineWidth set to 1.2
%         3. Default axes font set to Arial
%         4. Default text font set to Arial
%         5. Default text interpreter LaTeX (use $$ sign for math symbols)
%         6. Default legend interpreter LaTeX (use $$ sign for math symbols)
%         7. Default axes font size set to 10
%         8. Default line marker size set to 5
%         9. Default axes line width set to 1
%
% Copyright by Hewenxuan Li, 2021
% Last modified May 2, 2021
if nargin == 0
    options = 'Default';
end

if isequal(lower(options), 'default') 
    set(0,'DefaultFigureWindowStyle','docked')
    set(0, 'defaultFigureColor', 'w');
    set(0, 'defaultLineLineWidth', 1.2);
    set(0, 'defaultAxesFontName', 'Arial');
    set(0, 'defaultTextFontName', 'Arial');
    set(0, 'defaultTextInterpreter', 'tex');
    set(0, 'defaultLegendInterpreter', 'tex');
    set(0, 'defaultAxesFontSize', 10);
    set(0,'DefaultLineMarkerSize',5);
    set(0,'defaultAxesLineWidth',  1);

elseif isequal(lower(options), 'latex')
    set(0,'DefaultFigureWindowStyle','docked')
    set(0, 'defaultFigureColor', 'w');
    set(0, 'defaultLineLineWidth', 1.2);
    set(0, 'defaultAxesFontName', 'Arial');
    set(0, 'defaultTextFontName', 'Arial');
    set(0, 'defaultTextInterpreter', 'latex');
    set(0, 'defaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesFontSize', 10);
    set(0,'DefaultLineMarkerSize',5);
    set(0,'defaultAxesLineWidth',  1);

elseif isequal(lower(options),'elsevier')
    set(0,'DefaultFigureWindowStyle','docked')
    set(0, 'defaultFigureColor', 'w');
    set(0, 'defaultLineLineWidth', 1.2);
    set(0, 'defaultAxesFontName', 'Arial');
    set(0, 'defaultTextFontName', 'Arial');
    set(0, 'defaultTextInterpreter', 'tex');
    set(0, 'defaultLegendInterpreter', 'tex');
    set(0, 'defaultAxesFontSize', 10);
    set(0,'DefaultLineMarkerSize',5);
    set(0,'defaultAxesLineWidth',  0.5);
end
    defaultcolor = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980],	[0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
    mycolor = {'#0000b3','#e68a00','#e6b800','#e60000','#009933'};
    mycolor1 = {'#003399','#ff9900','#ffcc00','#ff3333','#00e64d'};
    blues = {'#000066','#000099','#0000cc','#0000ff','#3333ff','#6666ff','#9999ff'};
    pyblues = {'#cee0f2','#c2d9ee','#b1d2e8','#a0cbe2','#8bc0dd','#76b4d8','#62a8d3','#519ccc','#3282be','#2474b7','#1967ad','#0f59a3','#084c94','#083e80','#08306b'};
    matlabblues = {'#0071bc',' #007acc',' #008ae6',' #0099ff',' #1aa3ff',' #33adff',' #4db8ff',' #66c2ff',' #99d6ff',' #ccebff', ' #e6f5ff'};
    greys = {'#000000','#333333','#666666','#999999','#cccccc','#e6e6e6'};
    Mycolors.mycolor = mycolor;
    Mycolors.mycolor1 = mycolor1;
    Mycolors.blues = blues;
    Mycolors.greys = greys;
    Mycolors.defaultcolor = defaultcolor;
    Mycolors.pyblues = pyblues;
    Mycolors.matlabblues = matlabblues;
    Markers = {'O','+','x','*','s','d','v','^','>','<','p','h','O','+','x','*','s','d','v','^','>','<','p','h'};
    LineStyles = {'-','--','-.',':'};