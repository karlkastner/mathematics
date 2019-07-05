% Thu 27 Sep 16:47:27 CEST 2018
% Karl Kastner
%
%% dispersion determined by the hodges lehman method
%% asymptotic efficiency of dispersion estimates:
%% standard deviation:      E(s - hat s)/s = 2/sqrt(2 n) ~ 0.707/sqrt(n)			(100%)
%% hodges lehmann dispersion E(s-\hat s)/s = (pi/3)^2 /(sqrt(2 n)) ~ 0.775/sqrt(n)	(91%)
%% mad                      E(s-\hat s)/s ~ 1.17 s/sqrt(n)				(60%) 					
%% c.f. Shamos 1976
%% c.f. Bickel and Lehmann 1976
%% c.f. rousseeuw 1993
%% nb: rousseeuw uses the 25th percentile, which is more efficient for small sample sizes
function s = hodges_lehmann_dispersion(x)
	if (isvector(x))
		x = cvec(x);
	end

	s = zeros(1,size(x,2));
	for idx=1:size(x,2)
		XX = bsxfun(@minus,x(:,idx),x(:,idx).');
		XX = XX(triu(true(size(XX)),+1));
		s(idx)  = 1.0483*median(abs(XX));
	end
end

