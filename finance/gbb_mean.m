function gmu = gbb_geomean(S0,T,t,r,sigma);
	if (0)
	gmu = S0*exp(   log(gbm_mean(t,r,sigma,1)).*(1-t./T) ...
	              - log(gbm_mean(T-t,r,sigma,1)).*t./T ...
	         );
	end
	gmu = S0*ones(size(T));
%	B_ = exp(log(S).*(1-t./T) - t./T.*log(S(end,:)./(S.*S0)));
end

