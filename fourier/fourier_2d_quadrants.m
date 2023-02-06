function [f11,f12,f21,f22] = fourier_2d_quadrants(f)
	siz  = ceil(size(f)/2);
	f11  = f(1:siz(1),1:siz(2));
	f12  = f(1:siz(1),siz(2)+1:end);
	f21  = f(siz(1)+1:end,1:siz(2));
	f22  = f(siz(1)+1:end,siz(2)+1:end);
end

