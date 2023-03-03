% Fri  3 Mar 13:08:40 CET 2023
function [fc,Sc] = bandpass2d_pdf_mode(a,order,L,n)
%	Sl = lowpass2d_pdf_exact(fr,a,order);
%	Sb = 4*Sl.*(1.0-Sl);
%	Sl'*(1-Sl) - Sl*Sl' = 0
%	Sl' - 2 Sl*Sl' = 0 -> 1 - Sl = 0 min (1-Sl).^2
	fc0 = 1.0;
	fc = lsqnonlin(@(f) (bandpass2d_pdf(f,a,order)-1.0),fc0);
	if (nargout()>1)
	end
end

