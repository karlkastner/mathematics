% 2024-12-22 09:45:06.720724936 +0100
% Karl Kastner, Berlin
function [D1x,D1y] = fourier_derivative_coefficients(L,n)
	if (length(n)>1)
	[fx, fy] = fourier_axis_2d(L,n);
	else
		fx = fourier_axis(L,n);
		fy = 0; 		
	end
	ox = 2*pi*fx;
	oy = 2*pi*fy';
	D1x = 1i*ox;
	D1y = 1i*oy;
end

