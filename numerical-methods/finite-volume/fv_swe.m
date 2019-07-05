% Fri 17 Nov 09:20:12 CET 2017
%
%% wrapper for solving SWE
%
% % discretisation : number of points along domain
function [T X H U fv] = fv_swe(Ti,Xi,zbfun,wfun,cdfun,bc,icfun,Q0,a,opt)
	if (nargin()<10)
		opt = struct();
	end
	if (~isfield(opt,'limiter'))
		opt.limiter = 'upwind';
	end
	if (~isfield(opt,'dt_out'))
		% Tp         = 2*pi/omega;
		nt         = 1e3;
		opt.dt_out = (Ti(2)-Ti(1))/nt;
	end
	if (~isfield(opt,'QAflag'))
		opt.QAflag = true;
	end
	if (~isfield(opt,'nx'))
		opt.nx = 1e3;
	end
	switch (opt.limiter)
		case {'lax_friedrich'}
			fv = Lax_Friedrich();
		otherwise
			fv         = Reconstruct_Average_Evolve();
			fv.limiter = @(varargin) Flux_Limiter.(opt.limiter)(varargin{:});
	end % switch

	swe         = SWE();
	swe.QAflag  = opt.QAflag;

	% channel properties
	swe.wfun    = wfun;
	swe.zbfun   = zbfun;
	swe.cdfun   = cdfun;

	fv.pde      = swe;

	% initial condition
	% TODO the QAflag wrapper for the ic should be in SWE
	if (isstr(icfun))
		icfun = @(x) swe_ic(x, Xi, icfun, zbfun, wfun, cdfun, Q0, a);
	end
	fv.icfun       = @icfun_;

	% boundary conditions
	fv.bcfun       = bc;

	fv.dt_out      = opt.dt_out;
	fv.init(Xi,opt.nx);
	[T Y] = fv.solve(Ti);

	% extract
	X = fv.x;
	A = fv.pde.area(Y);
	H = fv.pde.depth(Y);
	U = fv.pde.velocity(Y);
	Q = fv.pde.discharge(Y);

	function y = icfun_(x)
		y = icfun(x);
		if (opt.QAflag)
			w = wfun(x);
			w = [w;w];
			y = w.*y;
		end
	end
end

