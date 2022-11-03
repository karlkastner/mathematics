% 2022-06-13 11:47:23.728711885 +0200
function y = bm_1d_fourier(n)
	e = randn(n(1),n(2));
	y = zeros(n(1),n(2));
	x = innerspace(0,1,n(1))';
	for idx=1:n(1)
		v = sqrt(2)*sin((idx-1/2)*pi*x);
		l = (idx-1/2)^2*pi^2;
		y = y + 1/sqrt(l)*v*e(idx,:);
	end
end
