% Thu 29 Jun 10:40:29 CEST 2017
% Karl Kastner, Berlin
%% transformation matrix for the short time fourier transform
% TODO set up by kronnecker product
function [A1, A2, DDl, DDr, A, obj] = stftmat(obj,time)
%	if (nargin()<6)
%		order = 1;
%	end
	if (max(obj.T) > obj.Ti)
		warning('Maximum period exceeds sample period');
	end

	obj.t0 = time(1);
	obj.tend = time(end);

	% number of samples
	ns = length(time);

	%t     = t0 + dt*(0:nt-1)';
	omega = 2*pi./obj.T;
	nf    = length(obj.T);

	% round to nearest integer multiple
	% TODO this is not good and should not come here

	%mt    = round(obj.Ti/obj.dt);
	%Ti    = mt*obj.dt;

	% number of unknowns per interval end point
	%ni = 2*nf+1;
	ni = obj.ni();

	% number of intervals
	%nti   = ceil(ns/mt)+1;
	nti  = obj.nti();
	%obj.nti = nti;

	% TODO, this should not be computed here
	%m_ = [ni, nti];

%	if (mt < ni)
%		% TODO only warn during fit, not during expansion
%		warning('More unknowns than samples');
%	end

	% matrix 1 : from frequency component at t to value at t
	% for each sample

	% row : index of sample point repeated for each component
	buf1 = flat(repmat((1:ns),ni,1));

	% colun : index of interpolated coefficient of sample point
	buf2 = flat(  repmat((0:ns-1)*ni,ni,1) ...
		     + repmat((1:ni)',1,ns) );
	buf3 = zeros(ns,ni);
	buf3(:,1) = 1;
	for idx=1:nf
		buf3(:,2*idx)   = cos(omega(idx)*time);
		buf3(:,2*idx+1) = sin(omega(idx)*time);
	end
	buf3 = flat(buf3.');

	% set up matrix
	A1   = sparse(buf1,buf2,buf3,ns,ns*ni);

	% interval index for each sample
	ii      = floor((cvec(time)-obj.t0)/obj.Ti)+1;

	% non-convex points
	fdx     = (ii<1) | (ii > nti);
	ii(fdx) = 1;

	% matrix 2 : from frequency component at support points to frequency components at time t
%	tic
	% sample point to frequency component per point sample point
%	buf1  = flat(repmat((0:ns-1)*ni,ni,1) + repmat((1:ni)',1,ns))';

	% frequency component is the same for all supporting points
%	buf1 = [buf1; buf1];
%	buf1 = flat(buf1);

%	buf1_ = flat(repmat((0:ns-1)*ni,1,1) + repmat(1,1,ns))';
%	buf1_ = flat([buf1_; buf1_]);

	% supporting coefficient
%	buf2 = flat(repmat((ii-1)'*ni,ni,1) + repmat((1:ni)',1,ns))';
%	buf2 = ([buf2; buf2+ni]);

%	norm(buf2-buf2_)
%	buf2 = flat(buf2);

%	buf2_ = flat(repmat((ii-1)',1,1) + repmat((1)',1,ns))';
%	buf2_ = flat([buf2_; buf2_+ni]);

	buf1 = flat(repmat((0:ns-1)*ni,ni,1) + repmat((1:ni)',1,ns))';
	buf2 = flat(repmat((ii-1)'*ni,ni,1) + repmat((1:ni)',1,ns))';

	switch (obj.order)
	case {-1} % do not use this, ill conditioned
		buf1 = repmat(buf1,2,1);
		buf2 = bsxfun(@plus,buf2_,ni*(0:1)');
		buf3 = 0.5*ones(size(buf2));
	case {0} % costant, y at 0
		buf3 = ones(size(buf2));
	case {1} % linear, y at -1,1 (for optimal conditioning)
		buf1 = repmat(buf1,2,1);
		buf2 = bsxfun(@plus,buf2,ni*(0:1)');

		% TODO use rem
		t      = (ii)' - rvec(time-obj.t0)/obj.Ti;
		t      = 1-2*t;
		buf31 = flat(repmat((1-t)/2,ni,1))';
		buf32 = flat(repmat((1+t)/2,ni,1))';
		buf3  = [buf31; buf32];
		buf3  = flat(buf3);

	case {3}
		% TODO, this fails, when there are less than 3 days

		% sample point index
		buf1 = repmat(buf1,4,1);

		% limit to interior points
		ii              = ii;
		ii(1 == ii)     = 2;
		ii(nti-1 == ii) = nti-2;

		% support point index
		buf2      = flat(repmat((ii-1)'*ni,ni,1) + repmat((1:ni)',1,ns))';
		buf2      = bsxfun(@plus,buf2,ni*(-1:2)');

		ti = obj.ti;
		id2 = (ii-1)';
		t = (rvec(time)-rvec(ti(id2)))./obj.Ti;
		
		buf3 =[];
		t = rvec(t);
		if (0)
			% lagrangian at 0..3 (not yet optimally conditioned)
			buf3(:,1) = flat(repmat(    - t.^3/6 + t.^2 - (11*t)/6 + 1,ni,1))';
			buf3(:,2) = flat(repmat(      t.^3/2 - (5*t.^2)/2 + 3*t,ni,1))';
			buf3(:,3) = flat(repmat(    - t.^3/2 + 2*t.^2 - (3*t)/2,ni,1))';
			buf3(:,4) = flat(repmat(          t.^3/6 - t.^2/2 + t/3,ni,1))';
		else % hermite y and y' at [-1.5 -0.5 0.5 1.5],
			% note, this is not optimally conditioned, optimally conditioned for t ~ 0.7381 t
			t = t-1.5;
			buf3(:,1) = flat(repmat( - t.^3/6 + t.^2/4 + t/24  - 1/16,ni,1))';
			buf3(:,2) = flat(repmat(   t.^3/2 - t.^2/4 - t*9/8 + 9/16,ni,1))';
			buf3(:,3) = flat(repmat( - t.^3/2 - t.^2/4 + t*9/8 + 9/16,ni,1))';
			buf3(:,4) = flat(repmat(   t.^3/6 + t.^2/4 - t/24  - 1/16,ni,1))';
%			buf3(:,1) = flat(repmat(     - t.^3/8 + t.^2/8 + t/8 - 1/8,ni,1))';
%			buf3(:,2) = flat(repmat(   t.^3/8 - t.^2/8 - (5*t)/8 + 5/8,ni,1))';
%			buf3(:,3) = flat(repmat( - t.^3/8 - t.^2/8 + (5*t)/8 + 5/8,ni,1))';
%			buf3(:,4) = flat(repmat(       t.^3/8 + t.^2/8 - t/8 - 1/8,ni,1))';
		end
		buf3 = buf3';
	otherwise
		error('here');
	end
	buf1 = flat(buf1);
	buf2 = flat(buf2);
	buf3 = flat(buf3);

	A2 = sparse(buf1,buf2,buf3,ns*ni,nti*ni);

	if (nargout() > 2)
		%D = derivative_matrix_1_1d(obj.nti,obj.tend-obj.t0,1);
		% the interval is normalised, so not real time
		Dr = derivative_matrix_1_1d(obj.nti,obj.nti-1,1);
		Dl = derivative_matrix_1_1d(obj.nti,obj.nti-1,-1);
		L = ones(obj.ni,1);
		if (0)
		for idx=1:length(obj.T)
			L(2*idx-1) = max(1,obj.Ti./obj.T(idx));
			L(2*idx)   = max(1,obj.Ti./obj.T(idx));
		end
		end
		DDr   = kron(Dr,diag(sparse(L)));
		DDl   = kron(Dl,diag(sparse(L)));
	end

	if (nargout() > 3)
		A = A1*A2;
	end

end % stftmat

function fA1()
	% allocate memory
	if (0)
	tic()
	buf1 = zeros(ns*ni,1);
	buf2 = zeros(ns*ni,1);
	buf3 = zeros(ns*ni,1);

	n = 0;
	for idx=1:ns
		% mean component
		n      = n+1;
		buf1(n) = idx;
		buf2(n) = (idx-1)*ni + 1;
		buf3(n) = 1;
		
		% frequency components
		for jdx=1:nf
			% cos
			n       = n+1;
			buf1(n) = idx;
			buf2(n) = (idx-1)*ni + 2*jdx;
			buf3(n) = cos(omega(jdx)*t(idx));
			% sin
			n       = n+1;
			buf1(n) = idx;
			buf2(n) = (idx-1)*ni + 2*jdx + 1;
			buf3(n) = sin(omega(jdx)*t(idx));
		end % for jdx
	end % for idx
	toc
	norm(buf1-buf1_)
	norm(buf2-buf2_)
	norm(buf3-buf3_)
	end % if 0

end

function fA2()
	if (0)
	% reallocation necessary, as A2 is larger
	tic
	no = 2;
	buf1 = zeros(ns*ni*no,1);
	buf2 = zeros(ns*ni*no,1);
	buf3 = zeros(ns*ni*no,1);
	n = 0;
	for idx=1:ns
		ldx = floor( (time(idx)-obj.t0)/obj.Ti ) + 1;
		%ldx = floor((idx-1)/mt) + 1;
		tl  = obj.t0+(ldx-1)*obj.Ti;
		tr  = tl+obj.Ti;
		% for each frequency component
		for jdx=1:ni
			n = n+1;
			buf1(n) = (idx-1)*ni + jdx;
			buf2(n) = (ldx-1)*ni + jdx;
			n = n+1;
			buf1(n) = (idx-1)*ni + jdx;
			buf2(n) = (ldx)*ni   + jdx;
			switch (obj.order)
			case {-1}
				buf3(n-1) = 1;
				buf3(n)   = 0;
			case {0} % constant interpolation
				% left
				buf3(n-1) = 0.5;
				% right
				buf3(n)   = 0.5;
			case {1} % linear interpolation
				% left
				buf3(n-1) = (tr - time(idx))/obj.Ti;
				% right
				buf3(n)   = (time(idx) - tl)/obj.Ti;
				%buf3(n-1:n)
				%pause
			otherwise
				error('here');
			end
		end % for jdx
	end
	toc
	norm(buf1-buf1_)
	norm(buf2-buf2_)
	norm(buf3-buf3_)
	%pause
	end

end
