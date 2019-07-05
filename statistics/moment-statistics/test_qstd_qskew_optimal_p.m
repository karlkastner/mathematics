% Mo 7. Mär 15:45:31 CET 2016
% Karl Kastner, Berlin
	n = [2e3,2e3];
	R = randn(n(1),n(2));

	m = 100;
	p = 0.5*(1:m)'/m;
	S = [];
	for idx=1:length(p)
		idx
		s  = qstd(R,1,p(idx));
		%s  = qskew(R,0,p(idx));
		%ms = mean(s)
		%s  = s/ms;
		S(idx,1) = std(s);
	end
	% reference is moment estimate
	%s0 = std(skewness(R));
	s0 = std(std(R));
	s_ = 1.4826*mad(R,1);
	mean(s_)
	s0(2) = std(s_);
	s_ = sqrt(2/pi)^-1*mad(R,0);
	mean(s_)
	s0(3) = std(s_);

	figure(1)
	clf
	plot(p,S)
	hline(s0)
	[Smin mdx] = min(S);
	p(mdx)
