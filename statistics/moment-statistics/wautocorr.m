% Wed  9 Aug 09:54:55 CEST 2017
% Karl Kastner, Berlin
%
%% autocorrelation for x if x is a vector, or indivvidually for the 
%% columns of x if x is a matrix
%% samples can be weighted
%%
%% c.f. box jenkins 2008 eq. 2.1.12
%%
%% c.f. autocorr_man4
%%
%% Note that it is faster to compute the acf in frequency space
%% as done in the matlab internal function
function [acorr avar] = wautocorr(w,x,k,type,resflag,biased)

	if (isvector(x))
		x = cvec(x);
	end
	n  = size(x,1);
	m  = size(x,2);

	if (nargin() < 3 || isempty(k))
		k = n-1;
	end
	if (nargin() < 4 || isempty(type))
		type = 'pearson';
	end
	if (nargin() < 5 || isempty(resflag))
		resflag = false;
	end
	if (nargin() < 6 || isempty(biased))
		% the default matlab definition is to use the biased estimate
		biased = true;
	end

	switch (type)
	case {'pearson'}
		if (~resflag)
			% residual with respect to the mean
			% TODO weighted mean
			mu = mean(x);
			d  = bsxfun(@minus,x,mu);
		else
			% the input is already a residual
			d = x;
		end
%		ww2 = sum(w.*w);
%		ww1 = sum(w);
%		w = bsxfun(@times,w,1./sum(w));
		% inner product with shifted data series
		sww0 = sum(w.^2);
		for idx=1:k
			ww = w(1:n-idx+1,:).*w(idx:end,:);
			if (~biased)
				%ww = bsxfun(@times,ww,1./sum(ww));
				avar(idx,:) = 1./sum(ww).*sum(ww.*d(1:n-idx+1,:).*d(idx:end,:));
				%avar(idx,:) = (n-idx+1)/n*sum(ww.*d(1:n-idx+1,:).*d(idx:end,:));
			else
				%w1  = w(1:n-idx+1,:);
				%w1  = bsxfun(@times,w1,1./sum(w1));
				%w2  = w(idx:end,:);
				%w2  = bsxfun(@times,w2,1./sum(w2));
				%w12 = w1.*w2;
				%w12 = bsxfun(@times,w12,1./sum(w12));
				%ww          = w(1:n-idx+1,:).*w(idx:end,:);
				%ww          = bsxfun(@times,ww,1./sum(ww));
				%avar(idx,:) = (n-idx+1)/n*sum(ww.*d(1:n-idx+1,:).*d(idx:end,:));
				avar(idx,:) = 1./sww0.*sum(ww.*d(1:n-idx+1,:).*d(idx:end,:));
			end
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
end % wautocorr

