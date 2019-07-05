% Do 16. Jul 09:49:49 CEST 2015
% Karl Kastner, Berlin
%
%% autocorrelation for x if x is a vector, or indivvidually for the 
%% columns of x if x is a matrix
%%
%% c.f. box jenkins 2008 eq. 2.1.12
%%
%% Note that it is faster to compute the acf in frequency space
%% as done in the matlab internal function
%
% TODO nan is problematic!
% function [acorr avar] = autocorr_man4(x,k,type,resflag,biased)
function [acorr avar] = autocorr_man4(x,k,type,resflag,biased)

	if (isvector(x))
		x = cvec(x);
	end
	n  = size(x,1);
	m  = size(x,2);

	if (nargin() < 2 || isempty(k))
		k = n-1;
	end
	if (nargin() < 3 || isempty(type))
		type = 'pearson';
	end
	if (nargin() < 4 || isempty(resflag))
		resflag = false;
	end
	if (nargin() < 5 || isempty(biased))
		% the default matlab definition is to use the biased estimate
		biased = true;
	end

	switch (type)
	case {'pearson'}
		if (~resflag)
			% residual with respect to the mean
			mu = mean(x);
			d  = bsxfun(@minus,x,mu);
		else
			% the input is already a residual
			d = x;
		end
		% inner product with shifted data series
		for idx=1:k
			avar(idx,:) = sum(d(1:n-idx+1,:).*d(idx:end,:));
			%N(idx,1) = size(d(1:n-idx+1,:),1);
			%fdx = isfinite(avar_);
			%l  = size(fdx,1);
			%l_ = sum(fdx);
			% if there are none-products, the sum has to be scaled by n/(n-n_an)
			% avar(idx,:) = (l./l_).*sum(avar_);

		end
		% normalisation
		if (biased)
			% biased, but lower mse, do not use for finite processes!
			avar = (1/n)*avar;
		else
			% unbiased, but higher mse
			%avar = bsxfun(@times,1./(n-1:-1:n-k)',avar);
			avar = bsxfun(@times,1./(n:-1:n-k+1)',avar);
		end

		% autocorrelation
		acorr = bsxfun(@times,avar,1./avar(1,:));
	otherwise
		error('not yet implemented');
%		a  = zeros(k,m);
%		for idx=1:k
%			for jdx=1:m
%				%ns = sum(d(1:n-idx+1,:) > d(idx:end,:));
%				%nl = sum(d(1:n-idx+1,:) < d(idx:end,:));
%				%a(idx,jdx) = corr(x(1:n-idx+1,jdx),x(idx:end,jdx));
%				a(idx,jdx) = corr(x(1:n-idx+1,jdx),x(idx:end,jdx),'type',type);
%			end
%		end	
	end
end % autocorr_man4

