% Wed 28 Jul 17:48:49 CEST 2021
% Karl Kästner, Berlin
%
%% standard deviation of a periodogram
%
function [sd] = periodogram_std(fx,S)
	fdx = (fx>=0);
	sd = wstd(S(fdx),fx(fdx));
end
