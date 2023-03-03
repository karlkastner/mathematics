% Fri 22 Apr 13:28:53 CEST 2022
% not normalized, max (S) = 1;
function Sb = bandpass2d_pdf(fr,a,order)
	% lowpass density
	Sl = lowpass2d_pdf(fr,a,1);
	% bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% higher order
	Sb = Sb.^order;
end

