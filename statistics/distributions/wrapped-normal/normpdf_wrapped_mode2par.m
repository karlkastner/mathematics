% Mon 13 Feb 13:42:27 CET 2023
function [fm,sf] = normpdf_wrapped_mode2par(fc,Sc)
	par0 = [fc,0.1./Sc];
	par = lsqnonlin(@resfun,par0);
	fm = par(1);
	sf = par(2);
	
function res = resfun(par)
	[fc_,Sc_] = normpdf_wrapped_mode(par(1),par(2));
	res = [fc_-fc,Sc_-Sc];
end

end

