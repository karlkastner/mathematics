% Thu 22 Jun 08:25:49 CEST 2017
%% initialize the interpolator with a set of sampling points
function [res, cflag, mse, obj] = init(obj,P0,val0,name,lambda,varargin)
	%if (isempty(obj.mesh))
		n  = round(length(P0)/10);
		%L  = range(P0);
		%xl = midrange(P0)-L/2;
		xlim = [min(P0),max(P0)];
		obj.remesh(xlim,n);	
	%end
	if (nargin()<4 || isempty(name))
		name = 'default';
	end
	if (nargin()<5 || isempty(lambda))
		lambda = obj.lambda;
	end
	
	% interpolate to point values
	[obj.vali.(name), res, cflag, mse, obj.msei.(name)] =  ...
	       obj.mesh.interp_tikhonov_1d(P0,val0,lambda,varargin{:});
end % init

