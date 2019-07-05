% Sat  9 Jul 15:26:11 CEST 2016
% Karl Kastner, Berlin
%
%% iverses wavelet transform for single frequency
%% (reconstruction of time series)
%% n : window lengths in multiples of filter period 1/f0
function [xs, x] = wavelet_reconstruct(w,dt,f0,n,winstr)

	if (nargin() < 5 || isempty(winstr))
		winstr = 'kaiser';
	end
%	if (nargin() < 6 || isempty(nanflag))
%		nanflag = false;
%	end
	if (isscalar(n))
		n = repmat(n,size(f0));
	end

	% undo phase correction
	t = dt*(0:length(w)-1)';
	%se = exp(2*pi*1i*(t)*f0);
	se = exp(2*pi*1i*(t+dt/2)*f0);
	w  = 1./se.*w;

	% allocate memory
	x = zeros(size(w));
	for idx=1:length(f0)
		% construct wavelet
		psi = wavelet(dt,f0(idx),n(idx),winstr);
%		psi = real(psi);

		% why has here real to be taken?
		wi = w(:,idx);
		wi = real(wi);

		% inverse transformation by convolution
		%xi   = conv([NaN;wi;NaN],psi,'same');
		xi   = conv_centered([NaN;wi;NaN],psi);

		% TODO, the admissibility constant is missing
		% scale, this is not required, as this lib does not apply forward scaling
		% x   = x.*f0.^2;

		% cut valid region
		x(:,idx) = xi(2:end-1);
	end % for idx
	% TODO there seems to be a bug in the matlab conv routine with flag 'same' set
	% this seems to be an odd-even problem and probably the fw and bw transform have each to be shifted by 1
	% see also + dt/2 fix in phase

%	x = circshift(x,[2 0]);
	x = circshift(x,[1 0]);
	x = real(x);
	xs = sum(x,2);
end % wavelet_reconstruct

