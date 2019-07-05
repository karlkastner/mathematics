% 2018-09-21 17:40:27.248294219 +0200
function y  = filt_hodges_lehman(x,n)
	m = floor((n+1)/2);
	y = zeros(size(x));
	n = length(x);
	for idx=1:n
		y(idx) = hodges_lehman(x(max(1,idx-m):min(n,idx+m)));
	end
end
