% syms x s; simplify(int(x.^2/s^2.*exp(-0.5*x.^2/s.^2)),'ignoreanalyticconstraints',true)
function a = rayleigh2d_area(s)
	a = pi*(2*pi)^(-3/2)*s;
end

