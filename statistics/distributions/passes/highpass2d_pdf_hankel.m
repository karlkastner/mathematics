% Fri 22 Apr 14:05:05 CEST 2022
% not normalized, lim_x->inf S(x) = 1;
function Sh = highpass2d_pdf_hankel(fr,a,order,varargin)
	Sl = lowpass2d_pdf_hankel(fr,a,[],varargin{});
	Sh = (1.0-Sl);
	if (nargin()>2 && ~isempty(order))
		Sh = Sh.^order;
	end
end

