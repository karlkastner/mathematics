% Mo 3. Aug 15:00:32 CEST 2015
% Karl Kastner, Berlin
%% approximate acf by the ar1 process
%% acf: autocovariance or autocorrerlation function
%% nf : skip first samples (for mixed geometric-arithmetic series (ARMA)
function [rho L delay] = ar1delay(acf,nf,finite)
	if (nargin() < 2)
		nf = 0;
	end
	if (nargin() < 3)
		finite = false;
	end
	rho = NaN;
	delay = NaN;
	n = length(acf);
	if (nf<n)
		delay=find(acf(nf+2:end) < acf(nf+1)*exp(-1),1,'first');
		if (~isempty(delay))
			% exact
			rho = (acf(nf+1+delay)/acf(nf+1))^(1/(delay+1));
		end
	else
		error('start must be smaller than series length');
	end
	if (finite)
		zlag = nf + find(acf(nf+1:end) < 0,1,'first');
%		rho  = fzero(@(rho) func(rho,n,zlag),rho);
%		rho
		%rho = lsqnonlin(@(rho) acf(nf+1:zlag-1)/acf(nf+1) - func2(rho,n,nf,zlag), rho);
		n = length(acf);
		W = (n:-1:1)';
		sW = sqrt(W);
		try
		rho = lsqnonlin(@(rho) sW(nf+1:end).*(acf(nf+1:end)/acf(nf+1) - func2(rho,n,nf,zlag)), 1/2*(1+rho));
		catch e
			e
			rho = NaN;
		end
%		rho
%		ratio0 = double(acf(nf+1+delay)/acf(nf+1))
%		rho = fzero(@(rho) ratio0 - ratio(rho,n,nf,delay),double(rho));
%		ratio(rho,n,nf,delay)
	end
	% Note, the effective sample size is n = 2*L (c.f. Leith)
	% 2L = n = -2/log(rho) ~ (1+r)/(1-r)
	L = -1./log(rho);
end

function func2 = func2(rho,n,nf,zlag)
	a     = acfar1(rho,n);
	func2 = a(nf+1:end)/a(nf+1);
%	func2 = a(nf+1:zlag-1)/a(nf+1);
end

function func = func(rho,n,delay)
	a = acfar1(rho,n);
	func = a(delay);
end

function ratio = ratio(rho,n,nf,delay)
	a = acfar1(rho,n);
	ratio = a(nf+1+delay)/a(nf+1);
end

