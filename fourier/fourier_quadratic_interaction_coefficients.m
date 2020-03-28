% Sun 15 Mar 14:06:14 +08 2020
% coefficients of the cubic interaction of a fourier series with itself
% 0,1,2,3,4
function a2 = fourier_quadratic_interaction_coefficients(a1,varargin)
	a2 = fourier_multiplicative_interaction_coefficients(a1,a1,varargin{:});
end

