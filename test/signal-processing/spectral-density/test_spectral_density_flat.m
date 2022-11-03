
% the spectral density of white noise is not defined,
% the magnitude of the periodogram goes with 2 L/n to zero
%
% the spectral density of the windowed white noise
% is white noise convolved with a gaussian filter, i.e. low passed white noise,
% the periodogram has also magnitude 2 L/n and the spectral density is also not defined

Lw=1; L = 100; n = 1e6; m=sqrt(n); x = linspace(-L/2,L/2,n)'; y = randn(n,1); w = normpdf(x,0,Lw); yy = [y,w.*y]; S = periodogram(yy,L); plot(meanfilt1(S,m)); hline(2*L/n); sum(yy.^2)*L/n

