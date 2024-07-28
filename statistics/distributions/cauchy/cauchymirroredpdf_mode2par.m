function [f0,s] = cauchywrappedpdf_mode2par(fc,Sc,par0)
	if (nargin()<3)
		[par(1),par(2)] = cauchypdf_mode2par(fc,Sc);
	end
	par = lsqnonlin(@(par) [fc,Sc] - fun(par),par);
	f0 = par(1);
	s  = par(2);
	function out = fun(in)
		[fc,Sc] = cauchywrappedpdf_mode(in(1),in(2));
		out = [fc,Sc];
	end
end
