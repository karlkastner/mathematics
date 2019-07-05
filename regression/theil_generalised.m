% 2015-11-23 11:56:57.763901349 +0100
% Karl Kastner, Berlin
%
%% generalization of the Theil-Senn operator to higher dimensions,
%% for arbitrary functions such as polynomials and multivariate regression
%% either higher order polynomials or glm
%% c.f. "On theil's fitting method", Pegoraro, 1991
%
function theil_multi(X, fit)
	% take m-tuples of points
	n = size(X,1);
	m = size(X,2);
	N = nchoosek(1:n,m);
	l = size(N,1);
	% sub2ind
	% N = bsxfun(@plus, N, l*(0:m-1));
	P = zeros(l,m);
	% fit parameter for each m-tupel
	parfor idx=1:l
		P(idx,:) = fit(X(N(idx,:),:));
	end
	% estimate parameter location
	p = median(P);
	% estimate parameter standard error
	sd   = median(abs(bsxfun(@minus,P,p)));
	% here contrary to Pegoraro, n is corrected down for the number of unknowns
	% TODO the error should better be estimated as some quantile difference of the parameter
	% note that the l-subsets are not independet
	serr = 1.86/sqrt(n-m)*sd;
end % theil_multi

