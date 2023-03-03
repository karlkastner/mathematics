% Fri 22 Apr 14:05:05 CEST 2022
% not normalized, lim_x->inf S(x) = 1;
function Sh = highpass2d_pdf(fr,a,order)
	Sl = lowpass2d_pdf(fr,a,order);
	Sh = (1.0-Sl);
end

