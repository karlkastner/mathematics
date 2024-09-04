% 2024-06-29 10:48:53.472138692 +0200
function R = laplaceacf(x,mu,s)
	R = cos(2*pi*x*mu)./(1 + s^2.*(2*pi).^2*x.^2);
end

