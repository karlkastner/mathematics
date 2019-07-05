% So 2. Aug 12:01:52 CEST 2015
%
%% pearson correlation coefficient
function rho = pearson(x,y,w)
	if (nargin() > 1 && ~isempty(y))
		mux = mean(x);
		muy = mean(y);
		dx  = x-mux;
		dy  = x-muy;
		if (nargin()<3)
			rho = sum(dx.*dy)./sqrt(sum(dx.*dx).*sum(dy.*dy));
		else
			% weighing
			rho = sum(w.*dx.*dy)./sqrt(sum(w.*dx.*dx).*sum(w.*dy.*dy));
		end
	else
		if (nargin() < 3 || isempty(w))
			C = cov_man(x);
		else
			C = cov_man(x,[],w);
		end
		Di = diag(sparse(1./sqrt(diag(C))));
		rho = Di*C*Di;
	end
end

