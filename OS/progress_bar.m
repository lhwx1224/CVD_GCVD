function progress_bar(i, tot, str);
%
%========================================================================================================
%========================================================================================================
%
% ==========    DESCRIPTION    ==========
% This function displays a sort of handmade progress bar on your command window.
%
% ==========    INPUT    ==========
% i:    (mandatory)     is the i-th iteration which give you an idea about the time progress (scalar);
% tot:  (mandatory)     is the last number to which i loop finishes (scalar);
% str:  (facultative)   this string is printed in command window to help you trace the calling function
%                       (string). If you use the mfilename command, you will see the string of the
%                       filename which is using the progress bar.
%
% ==========    EXAMPLES    ==========
% 1.\
% For instance if your loop is:
%       for j = 1:100
%           statement;
%       end
% you have to call the progress bar function as
%       progress_bar(j, 100)
% within you loop:
%       for j = 1:100
%           progress_bar(j, 100)
%           statement;
%       end
% 2.\
% You can implement the progress_bar function inside another function with a LOOP to solve:
%
%       function [varargout] = your_function_name(varargin);
%       ...
%       statements
%       ...
%       for i = 1:tot
%           progress_bar(i, tot, mfilename);
%           statements
%       end
%       statements
%
% With mfilename command you personalize the progress_bar function for your your_function_name function.
% You will see following stuff:
%
%                    your_function_name
% 0% #############################################_____ 100%
%                    94%
%
%========================================================================================================
%========================================================================================================

%-initialize
D = '#';
N = ' ';
p = {'0% '; ' 100%'};
par = 2:ceil(tot/50):tot;
if any(intersect(par,i))
    [nan c] = intersect(par,i);
    for z = 2:size(par,2)
        if z <= (size(par,2)/2-ceil(0.1*size(par,2)))
            N = [N ' '];
        end
        if z <= c
            D = [D '#'];
        elseif z > c
            D = [D '_'];
        end
    end
    %-print on command window
    clc
    if nargin > 2, disp([N str]); end
    disp([p{1} D p{2}]);
    disp([N num2str(ceil(i/tot*100)) '%']);
end
