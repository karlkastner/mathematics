function [a,b] = lowpass2d_Sckl2par(Sc,kl,par0)
	if (nargin<3)
	par0 = [1,1];
	end
	par = lsqnonlin(@(p) lowpass1dpdf([0,kl],p(1),p(2)) - Sc*[1,0.5],par0);
	a = par(1);
	b = par(2);
end

