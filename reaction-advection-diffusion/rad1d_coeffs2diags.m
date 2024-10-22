% Tue 17 Sep 14:58:51 CEST 2024
% Karl Kastner, Berlin
%% determine elements of the discretization matrix of a reaction-advection-diffusion
%% equation in one dimension
% TODO ax should be sampled at mid-points
function diags = rad1d_coeffs2diags(varargin)
	switch (length(varargin))
	case {2}
		ex = varargin{1}(:,1);
		ax = varargin{1}(:,2);
		re = varargin{1}(:,3);
		dx = varargin{2};
	case {6}
		ex = varargin{1};
		ax = varargin{2};
		re = varargin{3};
		dx = varargin{4};
	otherwise
		error('here');
	end
	% TODO diffusion:
	% e_+1/2 [y+1 - y0] - e_-1/2[y0 - y-1]
	% exmid = 0.5*(ex + right(ex,1))
	% axmid = 0.5*(ax + right(ax,1))

	% account for non-unit grid space
	ex = ex/(dx(1)*dx(1));
	ax = ax/dx(1);
	% sample ad mid-points
	% TODO, the coefficients should be directly passed at the grid points
	uex = 0.5*(ex+up(ex,1));
	dex = 0.5*(ex+down(ex,1));	
	% upwinding [0,1,-1], ax > 0
	%        or [1,-1,0]
	flagx = (ax>0);
	u         = uex; 
	u(~flagx) = u(~flagx) + ax(~flagx);
	c  = -(uex+dex) + re;
	c( flagx) = c( flagx) + ax( flagx);
	c(~flagx) = c(~flagx) - ax(~flagx);
	d         = dex;
	d(flagx)  = d(flagx)  - ax(flagx);

	diags     = [u(:),c(:),d(:)];
end

