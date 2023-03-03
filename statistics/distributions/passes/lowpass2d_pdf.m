% Fri 22 Apr 14:05:05 CEST 2022
% unnormalized (radial) density of the pth-order lowpass in two dimensions
% continuous space
function S = lowpass2d_pdf(fr,a,order)
	S = zeros(size(fr));
	for idx=1:numel(fr)
		S(idx) = integral(@(r) besselj(0,2*pi*r*fr(idx)).*r.*exp(-a*abs(r)),0,inf);
	end
	% S0 = integral(@(r) besselj(0,0).*r.*exp(-a*abs(r)),0,inf);
	S0 = 1./a^2;
	S = S/S0;
	if (nargin()>2)
		S = S.^order;
	end
end

