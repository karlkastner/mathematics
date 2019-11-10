% Sun 30 Sep 11:26:44 CEST 2018
function r2 = coefficient_of_determination(y,yp,np,mode)
	if (nargin() < 4 || isempty(mode))
		mode = 'pearson';
	end
	if (nargin() < 3 || isempty(np))
		np = 1;
	end
	switch (lower(mode))
		case {'pearson'}
			rho = corr(y,yp,'type','Pearson');
		case {'hodges-lehmann'}
			rho = hodges_lehmann_correlation(y,yp);
		case {'kendall'}
			rho = kendall_to_pearson(corr(y,yp,'type','Kendall'));
		case {'spearman'}
			rho = spearman_to_pearson(corr(y,yp,'type','Spearman'));
		otherwise
			error('here');
	end % switch
	
	% number of samples
	% TODO allow for effective sample size
	ns = length(y);
	% correct for sample size
	r2 = (ns-1)/(ns-np)*(1-rho^2);
end % coefficient_of_determination

