% So 2. Aug 11:59:21 CEST 2015
% Karl Kastner, Berlin
%
%% spearman's product moment coefficient
function rho = spearman(x,y,w)
	if (nargin() < 3)
		w = [];
	end
	[rx, tx]  = ranking(x,w);
	if (nargin()>1 && ~isempty(y))
		[ry ty]  = ranking(y,w);
	%	rho  = spearman_rank(rx,ry,tx,ty)
		rho = pearson(rx,ry);
	else
		rho = pearson(rx);
	end
end

