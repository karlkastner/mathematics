% TODO make elliptic for rectangular domains
function w= tukeywin_circular(n,pw)
	[fx,fy,frr] = fourier_axis_2d([1,1],n);
	f2 = n(1)/2;
	f1 = (1-pw)*f2;
	w = (1+cos(pi*(frr - f1)/(f2 - f1)))/2;
	%w = (1+cos(pi*(frr - pw*f2)/(f2 - pw*f2)))/2;
	w(frr < f1) = 1;
	w(frr > f2) = 0; % > or >=?
end
