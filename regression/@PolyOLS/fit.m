%Mon Nov 24 14:58:24 CET 2014
% Karl Kastner, Berlin
%% fit a polynomial function
%% like polyfit, but returns parameter error estimates
%% TODO automatically activate scaleflag
function [param, res, A, obj] = fit(obj,X,Y,W,ci)
	if (nargin() < 4)
		W = [];
	end
	if (isvector(W))
		W = cvec(W);
	end
	if (~isempty(W) && size(W,1) ~= size(Y,1))
		error('length of W must match length of Y');
	end

	% shift to avoind round off error (actually 0.5*(min+max))
	% better x-mid(minmax)/range
	if (obj.scaleflag)
		obj.x0 = mean(X);
		X      = X - obj.x0;
		obj.s  = max(X);
		X      = X / obj.s;
	else
		obj.x0 = 0;
		obj.s  = 1;
	end

	[param, A] = obj.fit_(X,Y,W,obj.order);

	% residual
	Yp     = A*param;
	res    = Yp-Y;
	np     = size(param,1);

	if (isempty(W)||size(W,3)>1)
		% sample size
		ne     = size(Y,1);
		% data variance
		var_y  = var(Y);
		% error estimate based on moments
		% variance of the residuals (aka mse)
		ssr    = sum(res.*res);
		serr2  = ssr/(ne-np);
	else
		if (isvector(W))
			W       = W/sum(sqrt(W)).^2;
			%ne     = sum(sqrt(W))/sum(W)

			% effective sample size
			%ne      = sum(sqrt(W)).^2/sum(W);
			% this is identical
			ne      = sum(W).^2/sum(W.^2); 
			n       = ne;
			ny      = size(Y,1);

			% var_y  = wvar(W_,Y);
			% variance
			W_      = repmat(W,1,size(Y,2));
			var_y   = sum(Y.*(W_.*Y)) - length(Y)*mean(sqrt(W_).*Y).^2;
			% mean square residual
			serr2   = sum(res.*(W_.*res));
			ssr     = ny*serr2;

			% sample size correction
			var_y   = (ny/(ne-1))*var_y;
			serr2   = (ny/(ne-np))*serr2;
			
			if (0)
			% error estimate based on moments
			% variance of the residuals (aka mse)
			W_     = repmat(W,1,size(Y,2));
			nw     = sum(W);
			ssr    = (ny/nw)*sum(res.*(W_.*res));
			serr2  = ssr/dof
			end
		else
			error('nondiagonal covmat not yet implemented');
		end
	end % if weighing matrix supplied
	
	% coefficient of correlation
	r2.pearson  = 1 - serr2./var_y;
	% identical to:
%	r2_ = corr(Y,Yp,'type','Pearson')^2;
%	r2.pearson_ = 1 - (ne-1)/(ne-np)*(1-r2_)

	if (obj.extended_statistics)
		%r2_ = spearman_to_pearson(corr(Y,Yp,'type','Spearman'))^2;
		%r2.spearman = (ne-1)/(ne-np)*(1-r2_);
		%r2_ = kendall_to_pearson(corr(Y,Yp,'type','Kendall'))^2;
		%r2.kendall  = (ne-1)/(ne-np)*(1-r2_);
		%r2_ = hodges_lehmann_correlation(Y,Yp);
		%r2.hodges_lehmann  = (ne-1)/(ne-np)*(1-r2_);
		obj.mad = median(abs(res));
	end

	% error covariance matrix
	% TODO compute from QR factors
	% TODO C needs to be weighed as well
	% TODO retransform parameters and errors
	C0     = inv(A'*A);
	if (1 == length(serr2))
		C = serr2*C0;
	else
		C = [];
	end

	if (obj.determine_leverage)
		obj.leverage = diag(A'*inv(A*A)*A);
	end

	obj.r2      = r2;
	

	obj.param   = param;
	obj.C       = C;
	obj.C0	    = C0;
	obj.params  = sqrt(diag(C0)*serr2);
	obj.nsample = ne;
	obj.serr    = sqrt(serr2);
	obj.ssr     = ssr;
end % fit

