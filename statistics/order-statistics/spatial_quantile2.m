% Sa 16. Jan 23:01:44 CET 2016
% Karl Kastner, Berli
%% spatial quantile
function Q = spatial_quantile(X,alpha,u)
%	u = u/norm(u);
%	u = (2*alpha-1)*u;
	% start value
	Q0 = mean(X,2);
	% has X to be shifted by median?
	opt.reltol = 1e-12; %sqrt(eps);
	opt.TolFun = 0; %sqrt(eps);
	opt.TolX = 1e-12; %sqrt(eps);
	Q = lsqnonlin(@(Q) Qfunc(Q,X,u), Q0,[],[],opt);
end

function Q = Qfunc(Q,X,u)
	X_min_Q = bsxfun(@minus,X,Q);
	Q = sum(Phi(u,X_min_Q));
end

function Phi = Phi(u,t)
	% for dimensions higher than 1, only 1-norm gives consistent quantiles
	p = 2;
	Phi = lpnorm(t,p) + (u'*t);
end

