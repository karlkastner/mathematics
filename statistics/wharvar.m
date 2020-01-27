% Fri  3 Jan 15:56:32 +08 2020
% variance of the weighted harmonic mean
% syms w x y; mu = wharmean([w,1-w],[x,y])
% dfdx = diff(mu,x)
% dfdx = w/x^2 mu^2
function s2 = wharvar(w,x)
	w = w/sum(w);
	mu = wharmean(w,x);
	s2 = ((w./x.^2).*mu.^2).^2.*(x-mu).^2;
	s2 = sum(s2);
end
