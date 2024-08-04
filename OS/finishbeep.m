function finishbeep(letter)
% finishbeep(letter) 
% 
% Syntax: finishbeep(letter)
% 
% Input: letter - string, a single string that is eigher an Latin
% alphabetic letter or a arabic numebr.
%
% Output: Morse Code Beeps of that letter
% -------------------- EXAMPLE -------------------------------------------
% String = 'A1 2';
%for i = 1:length(String)
%    finishbeep(String(i));
%    pause(0.2)
%end
%
% Hewenxuan Li 2021
% ------------------------------------------------------------------------
a = 1;
Tend = [0.5 0.5 1];
% Tend = 0.5;
fs = 1000;
freq = 261.6/2.5;
space = 0.22;
% ------------------------------------------------------------------------
%                               ALPHABET
% ------------------------------------------------------------------------
Dot = 0.5;
Dash = 1.2;

A = [Dot Dash];
B = [Dash Dot Dot Dot];
C = [Dash Dot Dash Dot];
D = [Dash Dot Dot];
E = [Dot];
F = [Dot Dot Dash Dot];
G = [Dash Dash Dot];
H = [Dot Dot Dot Dot];
I = [Dot Dot];
J = [Dot Dash Dash Dash];
K = [Dash Dot Dash];
L = [Dot Dash Dot Dot];
M = [Dash Dash];
N = [Dash Dot];
O = [Dash Dash Dash];
P = [Dot Dash Dash Dot];
R = [Dot Dash Dot];
S = [Dot Dot Dot];
T = [Dash];
U = [Dot Dot Dash];
V = [Dot Dot Dot Dash];
W = [Dot Dash Dash];
X = [Dash Dot Dot Dash]; 
Y = [Dash Dot Dash Dash];
Z = [Dash Dash Dot Dot];
% ------------------------------------------------------------------------
%                            ARABIC NUMBER
% ------------------------------------------------------------------------

One = [Dot Dash Dash Dash Dash];
Two = [Dot Dot Dash Dash Dash];
Thr = [Dot Dot Dot Dash Dash];
Four = [Dot Dot Dot Dot Dash];
Five = [Dot Dot Dot Dot Dot];
Six = [Dash Dot Dot Dot Dot];
Sev = [Dash Dash Dot Dot Dot];
Eit = [Dash Dash Dash Dot Dot];
Nin = [Dash Dash Dash Dash Dot];
Zero = [Dash Dash Dash Dash Dash];

if nargin == 0
    for i = 1:length(Tend)
        t = 0:1/fs:Tend(i);
        beep = a*sin(2*pi*freq*t);
        sound(beep)
        pause(space)
    end
else
    switch letter
        case 'A'
            Code = A;
        case 'B'
            Code = B;
        case 'C'
            Code = C;
        case 'D'
            Code = D;
        case 'E'
            Code = E;
        case 'F'
            Code = F;
        case 'G'
            Code = G;
        case 'H'
            Code = H;
        case 'I'
            Code = I;
        case 'J'
            Code = J;
        case 'K'
            Code = K;            
        case 'L'
            Code = L;
        case 'M'
            Code = M;
        case 'N'
            Code = N;
        case 'O'
            Code = O;
        case 'P'
            Code = P;
        case 'Q'
            Code = Q;
        case 'R'
            Code = R;
        case 'S'
            Code = S;
        case 'T'
            Code = T;
        case 'U'
            Code = U;
        case 'V'
            Code = V;
        case 'W'
            Code = W;
        case 'X'
            Code = X;
        case 'Y'
            Code = Y;
        case 'Z'
            Code = Z;
        case '1'
            Code = One;
        case '2'
            Code = Two;
        case '3'
            Code = Thr;
        case '4'
            Code = Four;
        case '5'
            Code = Five;
        case '6'
            Code = Six;
        case '7'
            Code = Sev;
        case '8'
            Code = Eit;
        case '9'
            Code = Nin;
        case '0'
            Code = Zero;
        case ' '
            Code = E;
            a = 0;
    end
    
    for i = 1:length(Code)
        t = 0:1/fs:Code(i);
        beep = a*sin(2*pi*freq*t);
        sound(beep)
        pause(space)
    end
end