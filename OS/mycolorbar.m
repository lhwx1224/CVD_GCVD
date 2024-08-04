function mycolormap = mycolorbar(MapName, Nscale)
% maps = mycolorbar(MapName, Nscale) renders the current graphics with customized
% colormap. When an output is designated, it outputs the colormap for
% future use.
% 
% Syntax: maps = mycolorbar(MapName, Nscale)
%
% Input: MapName - string, the name of the colormap, choose one from the
%                  list down below
%        Nscale - int, the total number of sample points in the color map,
%                 better if it is greater than 64.
% Output: maps - 2d array, RGB ratio triplet for the use of colormapping.
%
% MapNames: 'Viridis', 'BlackBody'
%
% Hewenxuan Li 2021
if nargin == 0
    MapName = 'Viridis';
    Nscale = [];
elseif nargin == 1
    Nscale = [];
end



switch lower(MapName)
    case 'viridis'
        % DEFINE THE VIRIDIS COLOR MAP
        color1 = '#FDE725';
        color2 = '#7AD151';
        color3 = '#22A884';
        color4 = '#2A788E';
        color5 = '#414487';
        color6 = '#440154';
        if isempty(Nscale)
            mycolormap = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6});
            if nargout == 0
                colormap(mycolormap);
            end
        else
            mycolormap = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6}, Nscale);
            if nargout == 0
                colormap(mycolormap);
            end
        end

    case 'blackbody'
        if isempty(Nscale)
            mycolormap = customcolormap(linspace(0,1,15), {'#000000', '#240f09', '#3e1611',...
                '#5a1b16', '#771e1a', '#96211e', '#b42622', '#c5411c', '#d65813', '#e47007',...
                '#e78d12', '#e9c327', '#e7de32', '#f6f090', '#ffffff'});
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap(linspace(0,1,15), {'#000000', '#240f09', '#3e1611',...
                '#5a1b16', '#771e1a', '#96211e', '#b42622', '#c5411c', '#d65813', '#e47007',...
                '#e78d12', '#e9c327', '#e7de32', '#f6f090', '#ffffff'}, Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
    case 'coolwarm'
        if isempty(Nscale)
            mycolormap = customcolormap(linspace(0,1,8), {'#3b4cc0','#6889ee','#9abaff', ...
                '#c9d8f0','#edd1c2','#f7a889','#e26a53','#b40426'});
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap(linspace(0,1,8), {'#3b4cc0','#6889ee','#9abaff', ...
                '#c9d8f0','#edd1c2','#f7a889','#e26a53','#b40426'}, Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
    case 'blues'
        if isempty(Nscale)
            mycolormap = customcolormap(linspace(0,1,8), {'#f7fbff','#dbe9f6','#bbd6eb', ...
                '#88bedc','#549ecd','#2a7aba','#0c56a0','#08306b'});
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap(linspace(0,1,8), {'#f7fbff','#dbe9f6','#bbd6eb', ...
                '#88bedc','#549ecd','#2a7aba','#0c56a0','#08306b'}, Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
    case 'rdbu'
        if isempty(Nscale)
            mycolormap = customcolormap(linspace(0,1,8), {'#67001f','#c1373a','#f09b7a', ...
                '#fbe3d5','#dceaf2','#87beda','#3079b6','#053061'});
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap(linspace(0,1,8), {'#67001f','#c1373a','#f09b7a', ...
                '#fbe3d5','#dceaf2','#87beda','#3079b6','#053061'}, Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
    case 'bentcoolwarm'
        if isempty(Nscale)
            mycolormap = customcolormap([0,0.142857143,0.285714286,0.428571429,0.571428571,0.714285714,0.857142857,1], ...
                fliplr({'#3b4cc0','#637dd5','#95ade3','#d1dcee','#ecd7cb','#de9e86','#cb634f','#b40426'}));
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap([0,0.142857143,0.285714286,0.428571429,0.571428571,0.714285714,0.857142857,1], ...
                fliplr({'#3b4cc0','#637dd5','#95ade3','#d1dcee','#ecd7cb','#de9e86','#cb634f','#b40426'}), Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
end