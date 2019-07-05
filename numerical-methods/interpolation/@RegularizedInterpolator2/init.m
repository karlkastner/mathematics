% Sun 23 Jul 12:22:47 CEST 2017
%% initialize the interpolator with a set of point samples
function [res cflag mse obj] = init(obj,P0,val0,name,lambda,varargin)
	if (isempty(obj.mesh))
		n  = round(length(P0)/100);
		xlim = [min(P0(:,1)),max(P0(:,1))];
		ylim = [min(P0(:,2)),max(P0(:,2))];
		obj.remesh(xlim,ylim,n,m);	
	end
	if (nargin()<4 || isempty(name))
		name = 'default';
	end
	if (nargin()<5 || isempty(lambda))
		lambda = obj.lambda;
	end
	
	% interpolate to point values
	if (obj.mseflag)
		[obj.vali.(name) res cflag mse obj.msei.(name)] =  ...
		       obj.mesh.interp_tikhonov_2d(P0,val0,lambda,varargin{:});
	else
		[obj.vali.(name) res cflag mse] =  ...
		       obj.mesh.interp_tikhonov_2d(P0,val0,lambda,varargin{:});
	end
end % init

