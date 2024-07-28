% for mirrored distribution
function R = normacf(x,fmu,fsd)
	t = 2*pi*x;
	R = cos(fmu*t).*exp(-0.5*fsd.^2*t.^2); 
end

