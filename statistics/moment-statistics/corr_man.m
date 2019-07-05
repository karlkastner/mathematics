% So 2. Aug 13:23:09 CEST 2015
% Karl Kastner, Berlin
%
%% correlation of two vectors
%
% TODO, at the moment kendall does not work for covariance matrices
function rho = corr_man(x,y,type,w)
	if (nargin() < 2)
		y = [];
	end
	if (nargin() < 3)
		w = [];
	end
	if (nargin()<3)
		type = 'pearson';
	end
	switch (lower(type))
	case {'pearson'}
		rho = pearson(x,y,w);
	case {'spearman'}
		rho = spearman(x,y,w);
	case {'kendall'}
		rho = kendall(x,y,w);
	otherwise
		error('unknown type for correlation coefficient');
	end
end % corr_man
