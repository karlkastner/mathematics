% Tue 19 Feb 20:36:56 CET 2019
%
% pdf of the skewed log distribution
% TODO mu is not mean in linear space, but log space, mgf?
% clf; x=innerspace(0.0,50,1e4); s=[-1,0,1]; d=[1]; for idx=1:length(s); for ddx=1:length(d); p = logskewpdf(x,1,d(ddx),s(idx)); plot(x,p); hold on; mu(idx)= sum(p.*x)/sum(x), end, end
function f = logskewpdf(x,mu,s,a)
	f = 2./(s*x).*normpdf((log(x)-mu)/s).*normcdf(a*(log(x)-mu)/s);
end

