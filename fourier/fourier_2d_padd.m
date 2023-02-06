function f = fourier_2d_padd(f,n)
	siz = size(f);
	% dissassemble
	[f11,f12,f21,f22] = fourier_2d_quadrants(f);
	% reassemble
	f = [f11, zeros(size(f11,1),n(2)-size(f,2)), f12;
	     zeros(n(1)-size(f,1),n(2));
	     f21, zeros(size(f21,1),n(2)-size(f,2)), f22];
end


