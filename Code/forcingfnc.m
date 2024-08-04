function [f, d, tspan, fs, dtn, sample] = forcingfnc(type, timing, finfo, rr)
% [F, D, TSPAN, FS, SAMPLE] = FORCINGFNC(TYPE,TIMING,FINFO,RR)
% 
% Syntax:
%       [f, d, tspan, fs, dt, sample] = forcingfnc(type, timing, finfo, rr)
%
% Input:
%       type - string, type of excitation (1) "Harmonic" or (2) "Random"
%       timing - structure, MATLAB structure who carries the timing
%       information of the forcing function; 
%           timing.Tstart - scalar, starting time of the forcing
%           timing.Tend - scalar, ending time of the forcing
%           timing.dT - scalar, time interval of the forcing (of the
%           function for interpolation);
%       finfo - structure, MATLAB structure which contains the forcing
%       information of the forcing function;
%           finfo.A - scalar, amplitude of the forcing
%           finfo.omega - scalar, frequency of the forcing 
%           finfo.windtype - string, type of window added to the forcing
%               function if the 'type' is "Random"
%           finfo.pwind - scalar (between 0 and 1), window size for the
%           window when the window type is "Harmonic" or "Rectangular"
%
% Output:
%       f - MATLAB structure, forcing described by spline
%       d - MATLAB strucutre, displacement described by spline
%       tspan - vector, time span of the generated forcing sample
%       (resampled if rr ~= 1) 
%       fs - scalar, sampling frequency after resampling
%       dt - scalar, sampling time after resampling 
%       sample - vector, sample of the designated forcing
%
% Examples:
% Ex1: Harmonic loading
%   f = A*cos(omega*t) where A = 0.3, omega = 1.3, no window; t0 = 0, t1 =
%   1000, dt = 0.2;
%   Input: timing.Tstart = 0; timing.Tend = 1000; timing.dT = 0.2;
%          finfo.A = 0.3; finfo.omega = 1.3; finfo.windtype = [];
%          finfo.pwind = [];
%          type = "Harmonic"; rr = 1 (no upsampling);
%   One calls the function:
%           [f, d, tspan, fs, sample] = forcingfnc(type, timing, finfo, rr)
%   
% Ex2: Random loading (Gaussian Random Process)
%   standard deviation = 0.3; t0 = 0, t1 = 1000, dt = 2; Rectangular
%   window; window effective range 90%; resampling rate = 100;
%   Input: 
%          Given fs
%          rr = 2;                         % Resampling Rate
%          type = "Random";                % Setup loading type 
%          timing.Tstart = 0;              % Starting time of simulation
%          timing.Tend = 1000;             % Ending time of simulation
%          timing.dT = 1/fs;               % Time interval (sampling time)
%          finfo.A = 0.3;                  % Standard deviation of the load
%          finfo.omega = [];               % Frequency of the load (empty)
%          finfo.windtype = "Rectangular"; % Window or not
%          finfo.pwind = 0.9;              % Window length
%          finfo.floc = "Base";            % Forcing Location
%          where fs/rr acts as an equivalent cut-off frequency
%   One calls the function:
%           [f, d, tspan, fs, sample] = forcingfnc(type, timing, finfo, rr)
%
% Programmed by Hewenxuan Li Oct, 2020
% Modified on Oct 22, 2020 - examples added
%             April 10, 2021 - added location of excitation and help.
%             April 29, 2021 - modified resampling machenism, focing
%             function will not be refitted, only the new time vector will
%             be returned. Help text updated.
%             May 1, 2021 - added burst random option in the random
%             section
% -------------------------------------------------------------------------

if nargin < 2
    timing.Tstart = 0;
    timing.Tend = 100;
    timing.dT = 0.1;
    finfo.A = 0.3;
    finfo.omega = 1.3;
    finfo.windtype = "Rectangular";
    finfo.pwind = 0.9;
    rr = 1;
elseif nargin < 3
    finfo.A = 0.3;
    finfo.omega = 1.3;
    finfo.windtype = "Rectangular";
    finfo.pwind = 0.9;
    rr = 1;
elseif nargin < 4
    rr = 1;
elseif nargin > 4 || nargin < 1
    error('Input number must be between 1 and 4!')
end

% Extract timing information
Tstart = timing.Tstart;
Tend = timing.Tend;
dt = timing.dT;
% Extract forcing information
A = finfo.A;
omega = finfo.omega;
windtype = finfo.windtype;
pwind = finfo.pwind;
floc = finfo.floc;
m = 1; % mass of the mass block
fs = 1/dt; % sampling frequency of the data for interpolation
t = Tstart:dt:Tend; % time vector
N = length(t); % number of samples in the data for interpolation
% ======================= Harmonic loading ================================
if isequal(type,"Harmonic")
    if isequal(floc,"base")
        sample = - m*A*omega^2*cos(omega*t); % Harmonic loading (in dimension of force)
    else
        sample = A*cos(omega*t); % Harmonic loading (in dimension of force)
    end
    d = spline(t,A*cos(omega*t)); % Displacement spline
    f = spline(t,sample); % Force spline
% ======================= Random loading ================================
elseif isequal(type,"Random")
%     rng('Default') % Set to default seed
    sample = m*A*randn(N,1); % Draw random sample from the Gaussian distribution
    sl = length(sample); % sample function size
    
    wind = zeros(size(sample));
    
    if isequal(windtype, "Harmonic")
        wind = sin(2*pi*pwind*t/Tend)'; % Create window
    elseif isequal(windtype, "Gaussian")
        wind = gausswin(sl);
    elseif isequal(windtype, "Rectangular")
        wind(floor(sl*pwind+1):floor(sl*(1-pwind))) = 1;
    elseif isequal(windtype, "BurstRandom")
        N = finfo.N;
        pwind = sl/N;
        for i = 1:2:N
            wind((i-1)*pwind+1:i*pwind) = 1;
        end
    else
        wind = ones(size(sample));
    end
    if size(wind,1)<size(wind,2)
        wind = wind';
    end
    sample = sample.*wind; % windowed sample function
    if isequal(floc,"base")
        d = spline(t,sample); % Displacement spline
        f = fnder(fnder(d)); % Force spline
    else
        f = spline(t,sample); % Force spline
        d = [];
    end
% ======================= Impulse loading ================================
elseif isequal(type,"Impulse")
    imp_width = pwind;        % window portion is used as pluse width here
    holdsample = 10;
    area = finfo.A;
    A = area/dt;
    Nimp = finfo.Nimp;
    imploc = fix(length(t)/Nimp);
    sample = 0*t;
    for i = 1:Nimp
    sample([(i-1)*imploc+(pwind-1) + holdsample,(i-1)*imploc+(pwind+1) + holdsample]) = A/imp_width;
    sample((i-1)*imploc+pwind + holdsample) = A/(imp_width/2);
    end
    if isequal(floc,"base")
        d = spline(t,sample);
        f = fnder(fnder(d)); % Force spline
    else
        f = spline(t,sample); % Force spline
        d = [];
    end
elseif isequal(type,"SI")
    rng('Default')
    % Harmonic loading
    sample = - m*A*omega^2*cos(omega*t); % Harmonic loading (in dimension of force)
%     d = spline(t,A*cos(omega*t)); % Displacement spline
%     f = spline(t,sample); % Force spline
    impact = zeros(size(sample));
    period = 2*pi/omega;
    n = floor(Tend/period); % number of full cycles
    p = floor(length(t)/n); % number of steps of each cycle
    Timp = floor(p*pwind);
    options = [];
    for i = 1+Timp:Timp:length(t)-Timp-1
%         impact(i) = 6*A*(rand(1,1)-0.5);
        impact(i:i+floor(Timp/2)) = A*(rand(1,floor(Timp/2)+1)-0.5);
    end
    f = spline(t, sample + impact);
    d = spline(t, A*cos(omega*t) + impact);
elseif isequal(lower(type), "freedecay")
    sample = 0*t; % No force applied
    f = spline(t,sample); % Force spline
    d = spline(t,sample); % Force spline
end
% New time: for ode or lsim function
% This part serves the purpose of a direct upsampling using the currently
% obtained forcing spline, e.g., if the the original time is,
%                     t = 0:0.01:10
% by imposing a resampling rate of 10, i.e., rr = 10; then, the new time
% is,
%                     tspan = 0:0.001:10
% Note that when the forcing is a sample from a random process, this rr is
% a resapling of the already fitted spline, not resampling  of the original
% random variable.
% If rr = 1, tspan returns the originally defined t vector.

dtn = dt/rr; % upsampling parameter
fs = 1/dtn; % upsampled sampling rate
tspan = Tstart:dtn:Tend; % upsampled time vector (will be returned)