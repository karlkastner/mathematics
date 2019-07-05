% Wed 26 Jun 10:26:56 CEST 2019
% Karl Kastner, Berlin
%
%% difference kernels for equispaced grids
%% c.f. Computing the Spectrum of the Confined Hydrogen Atom, Kastner, 2012
% 
% d     : order of the derivative (1,2 ...)
% order : order of accuracy
% 
% TODO implement higher order accuracy for other derivatives than 2
function [K] = difference_kernel(d,h,rep)
	if (nargin()<2 || isempty(h))
		h=1;
	end
	if (nargin()<3)
		rep = 0;
	end
	% first order derivative kernel
	K1  = 1/(2*h)*[-1, 0, 1];
	% second order derivative kernel
	K2  = 1./(h*h)*[ 1,-2, 1];
	switch (mod(d,2))
	case {1} % odd derivatives
		K = K1;
		for l=1:((d-1)/2)
			K = conv(K,K2);
		end		
	case {0} % even
		K = K2;
		% convolve K2 iteratively d/2 times, to get difference kernel
		% for derivative of order d
		for l=2:(d/2)
			K = conv(K,K2);
		end

		% higher order, only implemented for second order so far
		% note that this converges to the three-point difference
		% kernel for h->0
		if (2 == d)
			K2l   = K2;
			for l = 1:rep
				K2l = conv(K2l,K2);	
				c2l = nchoosek(2*l+1,l)*(l+1)^2;
				K   = [0,K,0] - (-1)^(l-1)/c2l*K2l;
			end
		end	
	end % switch mod(d,2)	
end % difference_kernel

