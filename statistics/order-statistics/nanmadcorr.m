% Wed  3 Aug 18:33:17 CEST 2016
% Karl Kastner, Berlin
%% proxy correlation by median absolute deviation
function [rn r] = madcorr(x,y)
	fdx = ~isnan(x) & ~isnan(y);
	x = x(fdx);
	y = y(fdx);

	mex = median(x);
	mey = median(y);
	madx = median(abs(x-mex));
	mady = median(abs(y-mey));
	r = 0.5*( median( abs((x-mex)/madx + (y-mey)/mady))  ...
		- median( abs((x-mex)/madx - (y-mey)/mady)) )
	% transformation under assumption of normality
	rn = r*sqrt(2-r^2);
	% limitation
	rn = max(-1,min(1,rn));
end

