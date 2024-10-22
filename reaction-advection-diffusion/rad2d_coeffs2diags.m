% Tue 17 Sep 14:58:51 CEST 2024
% Karl Kastner, Berlin
%
%% set up the elements of the discretization matrix of a two-dimensional
%% reaction-advection-diffusion equation, elements are stored in columns
%% corresponding to diagonals wrapped at the boundaries (circular boundaries)
% TODO pass coeffficients sampled at mid-points
function diags = rad2d_coeffs2diags(varargin)
	switch (length(varargin))
	case {2}
		% TODO, this needs some further thoughts how variables are stored
		ex = squeeze(varargin{1}(:,:,1,:));
		ey = squeeze(varargin{1}(:,:,2,:));
		ax = squeeze(varargin{1}(:,:,3,:));
		ay = squeeze(varargin{1}(:,:,4,:));
		% TODO combined reaction
		re = varargin{1}(:,5,:);
		dx = varargin{2};
	case {6}
		% n(1)*n(2)*nvar
		ex = varargin{1};
		ey = varargin{2};
		ax = varargin{3};
		ay = varargin{4};
		% n(1)*n(2)*nvar*nvar
		re = varargin{5};
		dx = varargin{6};
	otherwise
		error('here');
	end
	if (iscell(ex))
	% coupled equations
	for idx=1:length(ex)
		diags(:,:,idx) = rad2d_coeffs2diags(ex{idx},ey{idx},ax{idx},ay{idx},re{idx});
	end
	end
	else
		% single variable
		diags = rad2d_coeffs2diags(ex,ey,ax,ay,re,dx);
	end

function rad2d_coeffs2diags_(ex,ey,ax,ay,react,dx)

	% account for non-unit grid space
	ex = ex/(dx(1)*dx(1));
	ey = ey/(dx(2)*dx(2));
	ax = ax/dx(1);
	ay = ay/dx(2);

	nvar = size(ex,3);
%if (0)
%	ex = reshape(ex,n(1),n(2),nvar);
%	ey = reshape(ex,n(1),n(2),nvar);
%	ax = reshape(ex,n(1),n(2),nvar);
%	ay = reshape(ex,n(1),n(2),nvar);
%	% TODO 
%	re = reshape(re,n(1),n(2),nvar);
%end
	% TODO petrov-galerkin

	% sample ad mid-points
	% TODO, the coefficients should be directly sampled between grid-points
	lex = 0.5*(ex+left(ex,1));
	rex = 0.5*(ex+right(ex,1));	
	uey = 0.5*(ey+up(ey,1));
	dey = 0.5*(ey+down(ey,1));
	% upwinding [0,1,-1], ax > 0
	%        or [1,-1,0]
	flagx = (ax>0);
	flagy = (ay>0);
	u     = uey + ay.*(~flagy);
	%u(~flagy) = u(~flagy) + ay(~flagy);
	l  = lex + ax.*(~flagx);
	%l(~flagx) = l(~flagx) + ax(~flagx);
	%l(~flagx) = l(~flagx) + ax(~flagx);

	% reaction term
	c = react;
	%for idx=1:nvar
	% diagonal of diffusion terms (i=j)
	c(:,:,idx)  = -(lex+rex + uey+dey);
	%c( flagx) = c( flagx) + ax( flagx);
	%c(~flagx) = c(~flagx) - ax(~flagx);
	%c( flagy) = c( flagy) + ay( flagy);
	%c(~flagy) = c(~flagy) - ay(~flagy);
	% diagonal part of advection terms
	c(:,:,idx) = c(:,:,idx) + ax.*flagx - ax.*(~flagx);
	c(:,:,idx) = c(:,:,idx) + ay.*flagy - ay.*(~flagy);
	%end
	r         = rex - ax.*flagx;
	%r(flagx)  = r(flagx) - ax(flagx);
	d         = dey - ay.*flagy;
	%d(flagy)  = d(flagy) - ay(flagy);

	%diags     = [u(:),l(:),c(:),r(:),d(:)];
	diags     = [u(:),l(:),r(:),d(:),reshape(c,nn,[])];
	%if (nvar>1)
	%	diags = {u,l,r,d,c}
	%end
	end
end

