% Fri Feb 20 12:32:44 CET 2015
% Karl Kastner, Berlin
%% class for polynomial least squares
classdef PolyOLS < handle
	properties
		% error covariance matrix
		C
		C0
		% regression parameter
		param
		% standard error of regression parameters
		params
		% order of polynomial
		order 
		nsample
		x0;
		r2;
		serr;
		s;
		scaleflag;
		jk;
		ssr;
		leverage;
		determine_leverage;
		extended_statistics = true;
		mad;
	end % methods
	methods (Static) % prototypes
		slope     = slope(X,Y,W);
		[Y, slope] = detrend(X,Y,W);
		[param, A] = fit_(X,Y,W,order)
	end % methods (Static)
	methods
		function obj = PolyOLS(order,scaleflag)
			obj.order = order;
			if (nargin() > 1)
				obj.scaleflag = scaleflag;
			else
				scaleflag = false;
			end
		end
		function obj = jkfit(obj,varargin)
			obj.jk = Jackknife(@obj.regress_,obj.order);
			obj.jk.estimate(varargin{:});
		end
		function [Y S obj] = jkpredict(obj,X)
			func = @obj.predict_;
			[Y bias s2 C] = obj.jk.apply(func,X);
			S = sqrt(s2);
		end
%		function [varagout obj] = fit(obj,varargin)
%			[varargout obj] = obj.regress(varargin{:});
%		end
		function [nparam obj] = nparam(obj)
			if (isscalar(obj.order))
				nparam = obj.order+1;
			else
				nparam = sum(obj.order);
			end
		end
		% wrapper
		function rmse = rmse(obj,fdx)
			if (nargin()<2)
				rmse = obj.serr();
			else
				rmse = obj.serr(fdx);
			end
		end
		
	end % methods
end % classdef PolyOLS

