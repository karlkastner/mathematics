% Fri  8 Jul 15:08:33 CEST 2016
% Karl Kastner, Berlin
%
%% wavelet transform for single frequency
%% n : window lengths in multiples of filter period 1/f0
function w = wavelet_transform(x,dt,f0,n,winstr,pmin)
	x = cvec(x);
		
	if (nargin() < 5 || isempty(winstr))
		winstr = 'kaiser';
	end
	if (nargin() < 6 || isempty(pmin))
		pmin = 1;
	end

	w = zeros(length(x),length(f0));
	if (isscalar(n))
		n = repmat(n,size(f0));
	end

	if (pmin<1)
		fdx = isnan(x);
		x(fdx) = 0;
		o1 = double(~(fdx));
		o = double((fdx));
	end

	for idx=1:length(f0)
		% get wavelet
		[psi win sine] = wavelet(dt,f0(idx),n(idx),winstr);
	
		% transform by convolution
		if (1==pmin)
			wi   = conv_centered([NaN;x;NaN],conj(psi));
			% cut valid region
			w(:,idx)   = wi(2:end-1);
		else
			wi  = conv_centered(x,conj(psi));
			% weigh by number of valid samples
			s   = conv_centered(o1,win);
			s(s<pmin-sqrt(eps)) = NaN;
			% TODO the normalisation is bogus, this does not yield a valied wft transform with nans
			w(:,idx)  = wi./s;
		end
	end % for idx
	% rolling phase correction
	t = dt*(0:length(x)-1)';
%	se = exp(2*pi*1i*(t)*f0);
	se = exp(2*pi*1i*(t+dt/2)*f0);
	w  = se.*w;
end % wavelet_transform

