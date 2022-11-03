% 2021-09-28 17:46:58.462546616 +0200
% Karl Kästner, Berlin
%
%% qq-plot of a spectral density estimate by smoothing against the expected
%% beta-density
%
% function [q,p] = periodogram_qq(y,m)
function [q,p] = periodogram_qq(y,m)
	n = size(y,1);
	ni = round(n/m);
	% TODO use the mean-filtered periodogram here, as bartlett's density ratio
	% is not exactly betainv-distributed
	[S,Sd,SS] = periodogram_bartlett(y,1,m,ni);
 	SS = SS./S;
	ni_ = floor(ni/2);
	fdx = (1:ni_)';
	p=(1:m*ni_)'/(m*ni_+1);
 	q=[m*betainv(p,1,m-1),sort(flat(SS(fdx,:)))];
end

