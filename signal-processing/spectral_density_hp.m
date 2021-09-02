% Sat 26 Jun 21:04:27 CEST 2021
% spectral density of the high-pass
%function [S_h,S_h1] = spectral_density_hp(f,r0,fmax,order,varargin)
function [S_h,S_h1] = spectral_density_hp(f,r0,L,n,order,varargin)
	if (nargin()<4||isempty(order))
		order = 1;
	end
	[Sl,S_l1] = spectral_density_lp(f,r0,L,n,1,varargin{:});
	S_h = 1-S_l1.^order;
	% why this step?
	S_h = S_h./max(S_h);
	S_h = S_h.^2;
end 

