% Fri 22 Apr 13:28:53 CEST 2022
% not normalized, max (S) = 1;
function Sb = bandpass2d_pdf_hankel(L,n,a,order,varargin)
	% lowpass density
	Sl = lowpass2d_pdf_hankel(L,n,a,[],varargin{:});
	% bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% higher order
	if (nargin()>2 && ~isempty(order))
		Sb = Sb.^order;
	end
	% normalize
	df = 1/L;
	Sb = 2*Sb./(sum(Sb)*df);
end

