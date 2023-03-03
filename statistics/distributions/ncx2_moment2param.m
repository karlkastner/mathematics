% restricted parameter range!
% s^2 > 2 mu -> sd > sqrt(2*mu)
% there is also an upper limit
function [k,l] = ncx2_moment2param(mu,sd)
	% mu = k + l
	% s2 = 2k + 4l
	l = (sd.^2/2 - mu);
	k = mu - l;
end

