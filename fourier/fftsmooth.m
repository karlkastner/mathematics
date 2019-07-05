% Di 31. Mär 18:29:17 CEST 2015
% Karl Kastner, Berlin
%
%% smooth the fourier transform and determine upper and lower bound confidence intervals
%%
%% input :
%% f :
%% sfunc  : a smoothing function (for example fir convolution with rectangular window)
%%          returns filtered (mean) value and normalized fir window
%% nf     : window length
%% nsigma : number of standard deviations for confidnce intervals
%%
%%
%% output :
%% ff : filtered fourier transform
%% l : lower bound
%% u : upper bound
%
% TODO does this also apply to the periodogram?
% TODO shift f, so that low frequnecies are is in the centre
% TODO better pass window and directly convolve here
function [ff, fl, fu] = fftsmooth(f,sfunc,nf,nsigma)
	% filtered value and get filter window
	[ff, w] = sfunc(f,nf);
	% effective sample size (Kish)
	n = 2/sum(w.^2);
	c = normcdf(nsigma);
	% the fft bins are approximately independently chi^2 distributed
	fl = (n/chi2inv(c,n))*ff;
	fu = (n/chi2inv(1-c,n))*ff;
end

