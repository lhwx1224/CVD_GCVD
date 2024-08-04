function [poles, modes] = kITD(A)
% Ibrahim Time Domain Modal Identification
% [poles, modes] = ITD(Y) returns the discrete time poles (natural
% frequencies and dampling ratios) and the modes from a give set of
% responses (assumed to be free decay data sets).
% 
% syntax: [poles, modes] = ITD(Y)
%
% input: A - kernel for the ITD
%
% output: poles - 2d array, system poles in a matrix form
%         modes - 2d array, system modes in a matrix form
% 
% Number of records (response matrices from sensor output)
% Hewenxuan Li, Output-Only Modal Analysis Toolbox v0.0
% Created 2021 @ URI
% Last modified 09-20-2023 @ Cornell

% Check the data orientation, the default recognizes rows as different
% responses
if size(A,1) ~= size(A,2)
    error('Input A has to be square!')
end

[EVec, Eval] = eig(A);
poles = sqrt(Eval);
modes = EVec;