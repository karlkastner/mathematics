% 2021-09-28 17:46:58.462546616 +0200
% function [q,p] = periodogram_qq(y,m)
function [q,p] = periodogram_qq(y,m)
	n = size(y,1);
	ni = round(n/m);
	[S,Sd,SS] = periodogram_bartlett(y,1,m,ni);
 	SS = SS./S;
	ni_ = floor(ni/2);
	fdx = (1:ni_)';
	p=(1:m*ni_)'/(m*ni_+1);
 	q=[m*betainv(p,1,m-1),sort(flat(SS(fdx,:)))];
end

