% Fr 15. Jan 09:42:01 CET 2016
% Karl Kastner, Berlin
%
%% c.f. Oja 2008
%% is this the same as the oja simplex median (c.f. small 1990)?
function me = spatial_median(X)
	% damping, without damping no convergence
	% TODO if there are outliers, convergence is slow,
	% so better do a line search in direction of the gradient
	lambda = 1;
	abstol = sqrt(eps);
	n = size(X,2);
	m = size(X,1);
	% select mean as the start value
%	me = mean(X,2);
%	me = zeros(m,1);
	me = median(X,2);
	me_ = me;
	me__ = 0;
	% iterate
	while (1)	
		ss = zeros(m,1);
		no  = 0;
		X_min_me = bsxfun(@minus,X,me);
		ss = sum(spatial_sign(X_min_me),2);
		no = sum(sqrt(sum(X_min_me.*X_min_me)));
		me = me + lambda*(1/no)*ss;
		d = norm(me - me_)/lambda;
		d2 = norm(me - me__)/lambda;
		if (d2 < d)
			lambda = 0.5*lambda;
		end
		if (d < abstol)
			break;
		end
		me__ = me_;
		me_ = me;
	end
end

