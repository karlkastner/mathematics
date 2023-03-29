% Wed 27 Apr 10:55:29 CEST 2022
%
% lowpass filter the input x in the Frequency Domain
%
% TODO no need to provide dx, follows from size of x
function [y,S,R,r]=lowpass2d_ideal(x,L,dx,varargin)
	[S,R,r] = spectral_density_lowpass2d_ideal(size(x),dx,L,varargin{:});
	y = ifft2(sqrt(S).*fft2(x));
end

