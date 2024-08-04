function DataInspection(data, time, fs)
% DataInspection(data, time, fs) generates plot which contains (1) time
% series waveform and (2) the power spectral density function estimate
% using a Hanning window with fix length of 1/4 of the total length of the
% data.
% 
% input: data - data matrix 
%        time - time vector 
%        fs - sampling frequency
% 
% output: a figure with two subplots subplot (1): time waveform ; subplot
% (2): power spectral density estimate
% If one only input the data matrix, the plot generated will be in their
% sample-normalized frequency form.
%
% Hewenxuan Li, April 2021
figure
if size(data,1) < size(data,2)
    data = data';
end
windsize = round(size(data,1)/4);
if nargin == 3
    subplot(211)
    plot(time, data,'-')
    xlabel('Time (Seconds)')
    ylabel('Magnitude')
    title('Time Series')
    axis tight
    subplot(212)
    [pyy, f] = pwelch(data, hanning(windsize), windsize/2, [], fs);
    plot(f, 10*log10(pyy))
    xlabel('Frequency (Hz)')
    ylabel('Power/frequency (dB/Hz)')
elseif nargin == 2
    fs = 1/(time(2) - time(1));
    subplot(211)
    plot(time, data,'-')
    xlabel('Time (Seconds)')
    ylabel('Magnitude')
    title('Time Series')
    axis tight
    subplot(212)
    [pyy, f] = pwelch(data, hanning(windsize), windsize/2, [], fs);
    plot(f, 10*log10(pyy))
    xlabel('Frequency (Hz)')
    ylabel('Power/frequency (dB/Hz)')
else
    subplot(211)
    plot(data,'-')
    xlabel('Samples')
    ylabel('Magnitude')
    title('Time Series')
    axis tight
    subplot(212)
    pyy = pwelch(data, hanning(windsize), windsize/2);
    f = linspace(0, 1, size(pyy,1));
    plot(f, 10*log10(pyy))
    xlabel('Frequency (Hz)')
    ylabel('Power/frequency (dB/(rad/sample))')
end
end


