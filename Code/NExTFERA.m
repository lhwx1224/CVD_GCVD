% This function is combination of NExTF and ERA (NExTF is frequency-domain NExT and ERA is Eigensystem Realization Algorithim)
%Inputs :
%data: An array that contains response data.its dimensions are (nch,Ndata) where nch is the number of channels. Ndata is the total length of the data 
%refch: A vecor of reference channels .its dimensions (numref,1) where numref is number of reference channels
%window: window size to get spectral density
%N: Number of windows
%p: overlap ratio between windows. from 0 to 1
%fs: Sampling frequency 
%ncols: The number of columns in hankel matrix (more than 2/3*numref*(ceil(window/2+1)-1) )
%nrows: The number of rows in hankel matrix (more than 20 * number of modes)
%cut: cutoff value=2*no of modes
%shift: Shift value in the final row and column blocks (Increase EMAC sensitivity) usually =10
%EMAC_option: if this value equals to 1, EMAC will be independent of the number of columns (calculated only from observability matrix not from controllability) 
%Outputs :
%Result: A structure consist of the below components
%Parameters: NaFreq : Natural frequencies vector
%DampRatio: Damping ratios vector
%ModeShape: Mode shape matrix
%Indicators: MAmC : Modal Amplitude Coherence
%EMAC: Extended Modal Amplitude Coherence                                    
%MPC: Modal Phase Collinearity
%CMI: Consistent Mode Indicator
%partfac: Participation factor
%Matrices A,B,C: Discrete A,B and C matrices
function [Result] = NExTFERA(data,refch,window,N,p,fs,ncols,nrows,cut,shift,EMAC_option)
    IRF= NExTF(data,refch,window,N,p);                                %NExT
    inputs=length(refch);
    [Result]=ERA(IRF,fs,ncols,nrows,inputs,cut,shift,EMAC_option);    %ERA
end