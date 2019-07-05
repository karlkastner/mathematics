% Mi 2. Mär 08:47:29 CET 2016
% Karl Kastner, Berlin
%
% kernel density estimate in one dimension
% this is deprecated, as meanwhil kernel density estimates are available
% as standard functions in matlab
%
%% X   : ouput x axis bins
%% xi  : samples along x
%% m   : number of bins in X
%% fun : kernel function
%% pdf : propability density of xi
%
% function [pdf X] = kernel1d(X,xi,sx,m,fun)
function [pdf X] = kernel1d(X,xi,sx,m,fun)
	if (nargin() < 5 || isempty(fun))
		fun = @normpdf;
	end
%	if (nargin() < 2)
%		xi = linspace(min(X),max(X),round(sqrt(length(X))));
%	end
	if (isvector(xi))
		xi = cvec(xi);
	end
	n    = size(xi);
	if (isempty(X))
		if (nargin() < 4 || isempty(m))
			m   = round(sqrt(n(1)));
		end
		% create output matrices
		xmin = min(xi);
		xmax = max(xi);
		dx   = xmax-xmin;
		xmin = xmin-0.5*dx;
		xmax = xmax+0.5*dx;
		X(:,1) = linspace(xmin(1),xmax(1),m)';
	end
	if (isvector(X))
		X = cvec(X);
%		X = repmat(cvec(X),1,n(2));
	end
	m = size(X,1);

	
	% kernel density scale and covariance
	if (nargin() < 3 || isempty(sx))
		S = 8*std(xi)/sqrt(n(1));
	else
		S = sx;
	end
	if (size(S,1)~=n(1))
		S = repmat(S,n(1),1);
	end

	%Si = 1./S;

	% convolve samples with basis function
	% TODO this can indeed be much quicker implemented by convolution
	% TODO allow for wrapping around or mirroring at x
	pdf   = zeros(m,n(2));
	SS    = repmat(S,size(X,1),1);
	% expand basis function for each sample
	for jdx=1:n(2)
		for idx=1:n(1)
			pdf(:,jdx) = pdf(:,jdx) + fun(X,xi(idx,jdx),SS(idx,jdx));
		end
	end
	% normalise
	pdf = 1/n(1)*pdf;
	%end % for jdx

	if (nargout() < 1)
		plot(X,pdf);
	end
end % kernel1e

