% Mon 10 Aug 12:16:24 +08 2020
function [x_C,y_C,out] = bvp2c_stacked(odefun,bcfun,xi,nx,varargin)
	myexp = @(x) (1+x);
%	myexp = @exp;

	if (length(varargin)>0)
		inifun_C = varargin{1};
		opt = bvp2_check_arguments(varargin{2:end});
	else
		inifun_C = {};
		opt = bvp2_check_arguments(varargin{:});
	end

	% number of end-point coupled odes
	nc  = length(odefun);

	out = struct();

	% number of segments
	nxc = nx-1;
	
	% number of coupled odes
	cc  = feval(odefun{1},xi(1,1),1);
	neq = size(cc,3);

	nii = 0;
	nii_ = 0;
	mii = 0;
	for cdx=1:nc

		% segment end points
		x_C{cdx} = mesh1(xi(cdx,:),nx(cdx),opt.xs);

		% segment mid points
		xc_C{cdx} = mid(x_C{cdx});

		% segment lengths
		dx_C{cdx}  = diff(x_C{cdx});

		% number of equations, depending on order of ode is 1st (m=2) or 2nd order (m=3)
		mm(:,cdx)  = 2*ones(neq,1) + (0 ~= flat(cc(1,1,:)));
		% integrated number of equations
		mi(:,cdx) = cumsum([0;mm(:,cdx)]);
		%mii(cdx+1) = mii(cdx)+mi{idx};
		% number of points in each segment
		nii(cdx+1) = nii(cdx)+mi(end,cdx)*nxc(cdx);
		nii_(cdx+1) = nii_(cdx)+neq*nxc(cdx);
	end

	% allocate memory for differential operator
	AAA  = sparse([],[],[],nii(end),nii(end),4*nii(end));
	% allocate memory for inhomogeneous part
	bbb  = zeros(nii(end),1);

	% initial value of ypm
	% complex amplitude of the left and right going wave at segment mid points
	ypm    = zeros(nii(end),1);
	for cdx=1:nc
		if (~isempty(inifun_C))
			yc      = feval(inifun_C{cdx},xc_C{cdx});
			for idx = 1:neq
				% assign initial values into the inhomogeneous part
				% this does not matter, as yc is reassembled later
				% the inhomogeneous part is always the second part
disp('TODO');
% 				ypm(nxc*(mi(idx)+1)+1:nxc*(mi(idx)+2)) = yc((idx-1)*nxc+1:idx*nxc);
			end
		end
	end
	lll     = [];

	% solve non-linear system by picard iteration
	[ypm, out.cflag, out.kiter] = picard(@bvp2c_stacked_solve,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	y_C = inner2out_stacked(ypm);

% coupling condition : 
% mean flow :
%	Qup == Qdown
%	sign(x-xmid)*Q_i = 0

% first connecting channels :
%	sign(x-xmid)*Q_i = 0
% remaining connecting channels :
%	1/(i*o*w)*dQ/dx == 1/(i*o*w)*dQ/dx

% TODO, make general, this is customized for river tide computation
for cdx=1:length(ccon)
	ximid = 0.5*(xi(:,1)+xi(:,2));

	% channel ids
	cid = ccon(cdx).cid;
	% endpoint ids
	eid = ccon(cdx).eid;

	%
	p = ccon(cdx).eid;

	% for each parallel equation (frequency component)
	for idx=1:neq
		% row indices
	%	rid = ei(cdx,eid)+mi(idx) + TODO start/end
		% directions of channel with respect to bifurcation
		sig = sing(xi(cdx,eid) - ximid(cid));
	
		if (1 == idx)
			% mean flow
			% TODO
			% mean discharge for each channel, TODO, this has to go into bvp2c_assemble
			for cdx=1:nc
				% if Q0 is given, then set Q0, else
				% sum g hc^3 (z0_i - z0_i-1)/(xc_i - x_i-1) == - sum cd (Q0+Q1)^2_: rhs and unknown
				% when, w, h, cd constant, simplifies to g h^3/(w cd) (hn-h1)/(L-dx)
				for jd=
					A() = A() - g*hc./dxc;
					A() = A() + g*hc./dxc;
					% TODO, use friction coeffs
					A() =  cd./(wc.*hc.^2)*abs(Q0)
					b() = b() + 0.5*cd./(wc.*hc.^2).*abs(Q1).^2;
				end
			end
			% coupling condition
			A(rid(1),:) = 0;
			b(rid(1))   = 0;

			% eq 1 : sum Q0_i = 0
			for rdx=1:length(rid)
				A(rid(1),rid(rdx)) = 1;
			end
			b(rid(1)) = 0;

			% eq 2..nc : z0_i = z0_1
			for rdx=2:length(rid)
				A(rid(rdx),:)   = 0;
				A(rid(rdx)+1,:) = 0;
				b(rid(rdx),:)   = 0;
				b(rid(rdx)+1,:) = 0;
				A(rid(rdx),rid(1))       =  1;
				A(rid(rdx),rid(rdx))     = -1;
				A(rid(rdx)+1,rid(1)+1)   =  1;
				A(rid(rdx)+1,rid(rdx)+1) = -1;
			end
		else
			% tides
		% first row, equal discharge
		A(rid(1)+(1:3),:) = 0;
		b(rid(1)=(1:3))   = 0; % no external inflow/outflow at bi
		for rdx=1:length(rid)
			% f consits of 3 elements (columns) and is shifted
			% Q-, Q+ and constant part
			A(rid(rdx),rid(rdx))     = sig(rdx);
			A(rid(rdx)+1,rid(rdx)+1) = sig(rdx);
			A(rid(rdx)+2,rid(rdx)+2) = sig(rdx);
		end
		% reamining, : z_i == z_1
		for rdx=2:length(rid)
			A(rid(rdx),:)   = 0;
			A(rid(rdx)+1,:) = 0;
			A(rid(rdx)+2,:) = 0;
			b(rid(rdx))     = 0;
			b(rid(rdx)+2)   = 0;
			b(rid(rdx)+3)   = 0;
			% z- = 1/(iow) dQ/dx = 1/(iow) l Q
		%	A(rid(rdx),rid(1))     =  dQm/dx TODO
		%	A(rid(rdx),rid(rdx)) = -dQm/dx TODO
			% TODO, what about constant?
			% z+
		%	A(rid(rdx)+2,rid(1)+2) =  dQp/dx TODO
		%	A(rid(rdx)+2,rid(1)+2) = -dQp/dx TODO
		end
	end
end

-> coupling for z0
-> coupling bc of z0_i == z0_j (n-1 conditions)
		  sum Q = 0	(nth-condition)
-> equations for unknown Q0 (nxc+1) :
-> replace this row by inflow, where Q0 is given
-> watch out in iteration, initial value of Q0 has to be non-zero



function ypm = bvp2c_stacked_solve(ypm)	
	nci = cumsum([0,rvec(nxc)]);
	for cdx=1:nc
		[AA,rr,ll] = bvp2c_assemble(    ypm(nii(cdx)+1:nii(cdx+1)), ...
						odefun{cdx}, ...
					        bcfun{cdx}, ...
						xi(cdx,:), ...
						xc_C{cdx}, ...
						dx_C{cdx}, ...
						neq, ...
						nxc(cdx), ...
						mm(:,cdx), ...
						mi(:,cdx) ...
					   );
		AAA(nii(cdx)+1:nii(cdx+1),nii(cdx)+1:nii(cdx+1)) = AA;
		rrr(nii(cdx)+1:nii(cdx+1),1) = rr;
		lll(nci(cdx)+1:nci(cdx+1),1:2,:) = ll;
	end
	% solve
	ypm = AAA \ rrr;
end

	function y_C = inner2out_stacked(ypm)
		%y = zeros(sum(nxc)*neq,1);
		nxi = cumsum([0,rvec(nx)]);
		nci = cumsum([0,rvec(nxc)]);
		for cdx=1:nc
		   for id=1:neq
		     r = lll(nci(cdx)+1:nci(cdx+1),1,id);
		     if (2 == mm(id,cdx))
			if (0 ~= r(1,1,id))
				y_ = (   ypm(nii(cdx)+(mi(id,cdx)*nxc(cdx)+1:mm(id,cdx):mi(id+1,cdx)*nxc(cdx))) ...
				       + ypm(nii(cdx)+(mi(id,cdx)*nxc(cdx)+2:mm(id,cdx):mi(id+1,cdx)*nxc(cdx))) );
			else
			% degenerated linear function
			y_ = ypm(mi(id)*nxc+2:mm(id):mi(id+1)*nxc);
			end
			%y(nx*(id-1)+1:nx*id) = inner2outer(y_);
			y_C{cdx}(:,id) = inner2outer(y_);
		     else
			y_ = (   ypm(nii(cdx)+(mi(id)*nxc(cdx)+1:mm(id,cdx):mi(id+1,cdx)*nxc(cdx))) ...
			       + ypm(nii(cdx)+(mi(id)*nxc(cdx)+2:mm(id,cdx):mi(id+1,cdx)*nxc(cdx))) ...
			       + ypm(nii(cdx)+(mi(id)*nxc(cdx)+3:mm(id,cdx):mi(id+1,cdx)*nxc(cdx))) ...
		             );
			y_C{cdx}(:,id) = inner2outer(y_);
			%y(nxi(cdx)+(nx(cdx)*(id-1)+1:nx(cdx)*id)) = inner2outer(y_);
		     end
		   end % for id
		end % for cdx
	end % inner2outer

end % function bvp2c_stacked

