% Fri 22 Apr 13:28:53 CEST 2022
% not normalized, max (S) = 1;
function Sb = bandpass2d_pdf_hankel(fr,a,order,varargin)
	% lowpass density
	Sl = lowpass2d_pdf_hankel(fr,a,[],varargin{:});
	% bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% higher order
	if (nargin()>2 && ~isempty(order))
		Sb = Sb.^order;
	end
	% TODO normalize
end

