% Do 16. Jul 09:49:49 CEST 2015
% Karl Kastner, Berlin
%
%% compute autocorrelation for x if x is a vector, or indivvidually for the 
%% columns of x if x is a matrix
%
%% box jenkins 2008 eq. 2.1.12
%% TODO nan is problematic!
%
%% Note that it is faster to compute the acf in frequency space
%% as done in the matlab internal function
function [a dd] = autocorr_man4(x,k,mode,resflag)

	if (isvector(x))
		x = cvec(x);
	end
	n  = size(x,1);
	m  = size(x,2);

	if (nargin()<2 || isempty(k))
		k = n-1;
	end
	if (nargin() < 3 || isempty(mode))
		mode = 'pearson';
	end
	if (nargin() < 4 || isempty(resflag))
		resflag = false;
	end

	if (strcmp(lower(mode),'pearson'))
		if (~resflag)
			% residual with respect to the mean
			mu = nanmean(x);
			d  = bsxfun(@minus,x,mu);
		else
			% the input is already a residual
			d = x;
		end
		% inner product with shifted data series
		% TODO this can be computed by nanconv
		for idx=1:k
			dd_ = d(1:n-idx+1,:).*d(idx:end,:);
			fdx = isfinite(dd_);
			l  = size(fdx,1);
			l_ = sum(fdx);
			% if there are none-products, the sum has to be scaled by n/(n-n_an)
			dd(idx,:) = (l./l_).*nansum(dd_);
		end
		% normalise
		% note: normalisation with n is actually not necessary,
		% as it drops out during normalisation with a1
		dd = dd/n;
		a = bsxfun(@times,dd,1./dd(1,:));
	else
		a  = zeros(k,m);
		for idx=1:k
			for jdx=1:m
				%ns = sum(d(1:n-idx+1,:) > d(idx:end,:));
				%nl = sum(d(1:n-idx+1,:) < d(idx:end,:));
				%a(idx,jdx) = corr(x(1:n-idx+1,jdx),x(idx:end,jdx));
				a(idx,jdx) = corr(x(1:n-idx+1,jdx),x(idx:end,jdx),'type',mode);
			end
		end	
	end
end % acf_man

