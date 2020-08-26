% Sat 28 Oct 14:43:06 CEST 2017
%
% function [AAA,bbb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni,nci,npi,npii,jfun,dischargeisvariable)
%
%
% output :
%	AAA : left hand side discretization matrix
%	bbb : right hand side vector
%	out(id).lll : eingevalues of the characteristic functions
%
% TODO, treatment of degenerated second order equations a y'' + by' + cy = 0, with c = 0
% TODO, treatment of degenerated first order solution at junctions
function [AAA,bbb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni,nci,npi,npii,jfun,dischargeisvariable,bcarg)
	nc = length(odefun);

	% allocate memory for differential operator
	% TODO global buffers
	AAA  = sparse([],[],[],npii(end,end)-1,npii(end,end)-1,6*npii(end,end)-1);

	% allocate memory for inhomogeneous part
	bbb  = zeros(npii(end,end)-1,1);

	% for each edge in graph (channel in network)
	for cdx_=1:nc
		% set up discretisation matrix for coupled odes along edge (channel)
		[AA,bb,ll,cc] = bvp2c_assemble_AA( ...
                                                  ypm(npii(1,cdx_):npii(end,cdx_)-1) ...
						, odefun{cdx_} ...
					        , bcfun{cdx_} ...
						, xi(cdx_,:) ...
						, out(cdx_).xc ...
						, out(cdx_).dx ...
						, neq ...
						, nxc(cdx_) ...
						, oo ...
						, ni(:,cdx_) ...
						, nci(:,cdx_) ...
						, npi(:,cdx_) ...
					   );
		% write to global discretization matrix
		AAA(npii(1,cdx_):npii(end,cdx_)-1,npii(1,cdx_):npii(end,cdx_)-1) = AA;
		bbb(npii(1,cdx_):npii(end,cdx_)-1,1)     = bb;
		out(cdx_).ll = ll;
		out(cdx_).cc = cc;
		%.lll(npii(1,cdx_):npii(end,cdx_)-1,1:2,:) = ll;
	end % for cdx_

	if (nargin()>12 && ~isempty(jfun))
		couple_junctions();
	end

function couple_junctions()
    if (~dischargeisvariable)
		error('here')
    end
    % mid point of channels (for determining direction)
    ximid = 0.5*(xi(:,1)+xi(:,2));

    % TODO, make general, this is semi-customized for river tide computation
    % for each junction
    for jdx=1:length(jfun)

	% channel direction
	% dir = sign(xi(:,2)-xi(:,1));
    
    	% cid : channel indices
	% eid : endpoint indices {1,2},
	% p   : scales, for matching level of the frequency compontents,
	%               these are 1/(iow) 
    	[cid, eid, p] = feval(jfun{jdx});

	% fetch end point properties
	% segment length
	dx = zeros(length(cid),1);

	% eigenvalues
	l  = zeros(length(cid),2,neq);
	% indices
	rid_ = zeros(neq,length(cid));
	col_ = zeros(neq,length(cid));
	for cdx=1:length(cid)
		switch (eid(cdx))
		case {1}
			dx(cdx)     = out(cid(cdx)).dx(1);
			l(cdx,:,:)  = out(cid(cdx)).ll(1,:,:);
			rid_(:,cdx) = npii(1:end-1,cid(cdx));
			col_(:,cdx) = npii(1:end-1,cid(cdx));
		case {2}
			dx(cdx)     = out(cid(cdx)).dx(end);
			l(cdx,:,:)  = out(cid(cdx)).ll(end,:,:);
			rid_(:,cdx) = npii(2:end,cid(cdx))-1;  
			% TODO this is a quick and dirty hack
			% nb : tidal shift by 3, because there are 3 parts
			%      mean shift by 3, because there are 2 parts and the mean discharge
			col_(:,cdx) = npii(2:end,cid(cdx)) - 3.*(1==cvec(oo)) - 3.*(2==cvec(oo));
			%rid_(:,cdx) = rid_(:,cdx) - 1;
		otherwise
			error('here');
		end
    		dir(cdx) = sign(xi(cdx,eid(cdx)) - ximid(cid(cdx)));
	end % for cdx

    	% for each parallel equation (frequency component)
    	for edx=1:neq
    		% row indices
		rid = rid_(edx,:)';
		col = col_(edx,:)';

		switch (oo(edx))
		case {1} % match mean frequency component
			% condition for first connecting channel :
    			% conservation of (tidally averaged) discharge
			% sum dir*Q0_i  == 0 (1 equation)
			% reset row
			AAA(rid(1),:) = 0;
			bbb(rid(1))   = 0;
			for cdx=1:length(cid)
				% mean discharge Q0 has only 1 variable
				%AAA(rid(1),rid(cdx)) = dir(cdx);
				% here, the column indices are for the tidally averaged discharge,
				% not of the tidally averaged water level (rid)
				% TODO avoid magic number for index "2"
				AAA(rid(1),npii(2,cid(cdx))-1) = dir(cdx);
			end % for cdx
				
			% condition for remaining channels meeting at junction:
			% continuty of surface elevation (equal tidally averaged water level)
    			% z0_i - z0_1 == 0, i>1 (n-1 equations)
			% the water level does not depend on channel direction,
			% so dir is not factored in
    			for cdx=2:length(cid)
				% reset row
				AAA(rid(cdx),:) = 0;
				bbb(rid(cdx))   = 0;
			% note : solution for water level is always degenerated
				% homogeneous part	
% TODO check cc (!)
				if (0) %abs(l(1,1,edx)) ~= 0)
					% this actually never happens
					AAA(rid(cdx),col(1))     = -myexp(dir(1)*0.5*l(1,1,edx).*dx(1));
				else
					AAA(rid(cdx),col(1))     = -dir(1)*0.5*dx(1);
				end
				% inhomogeneous part
				AAA(rid(cdx),col(1)+1)   = -1;

				% homogeneous part	
% TODO check cc (!)
				if (0) %abs(l(cdx,1,edx)) ~= 0)
					% this actually never happens
					AAA(rid(cdx),col(cdx))   = +myexp(dir(cdx)*0.5*l(cdx,1,edx).*dx(cdx));
				else
					AAA(rid(cdx),col(cdx))   = +dir(cdx)*0.5*dx(cdx);
				end
				% inhomogeneous part
				AAA(rid(cdx),col(cdx)+1) = +1;
    			end % for cdx
		case {-1}
			% TODO deprecated
			% nothing to do (Q0 is matched for 1==edx)
		case {2}
			% frequency components (tide)

			% condition for first connecting channel :
			% conservation of tidal discharge of e-th frequency component
			% sum s Qt_i == 0, i = 1 .. n

			% reset row
    			AAA(rid(1),:) = 0;
			% no external inflow/outflow at junction
    			bbb(rid(1))   = 0;
    
    			for cdx=1:length(rid)
    				% discharge consits of 3 unknowns (columns)
				% left going part (Q-)
    				AAA(rid(1),col(cdx))   = dir(cdx)*myexp(dir(cdx)*0.5*l(cdx,1,edx)*dx(cdx));
				% inhomogeneous part
    				AAA(rid(1),col(cdx)+1) = dir(cdx);
				% right going part (Q+)
				AAA(rid(1),col(cdx)+2) = dir(cdx)*myexp(dir(cdx)*0.5*l(cdx,2,edx)*dx(cdx));
    			end % for cdx

			% condition for remaining connecting channels :
			% equal (oscillation) of water level of e-th frequency component
			%    zt_i - zt_1 == 0
			% => 1/(iow_i) dQt_i/dx - 1/(iow_1) dQ_1/dx == 0
			% => 1/(iow) (l_i^- Q_i^- + l_i^+ Q_i^+) - 1/(iow_1) (l_i^- Q_i^- + l_i^+ Q_i^+) == 0
			% p = 1/iow, 1/io can be factored out, but 1/w not
			% sign of channel direction is not factored in, as the water level is not a directional vector quantity
			for cdx=2:length(rid)
				% reset row
    				AAA(rid(cdx),:) = 0;
    				bbb(rid(cdx))   = 0;

				% left going
    				AAA(rid(cdx),col(1))     = -p(1)*l(1,1,edx)*myexp(dir(1)*0.5*l(1,1,edx)*dx(1));
				% derivative of inhomogeneous part is zero
				% right going
				AAA(rid(cdx),col(1)+2)   = -p(1)*l(1,2,edx)*myexp(dir(1)*0.5*l(1,2,edx)*dx(1));

				% left going
    				AAA(rid(cdx),col(cdx))   = +p(cdx)*l(cdx,1,edx)*myexp(dir(cdx)*0.5*l(cdx,1,edx)*dx(cdx));
				% derivative of inhomogeneous part is zero
				% right going
				AAA(rid(cdx),col(cdx)+2) = +p(cdx)*l(cdx,2,edx)*myexp(dir(cdx)*0.5*l(cdx,2,edx)*dx(cdx));
			end % for cdx
		otherwise
			error('here');
		end % switch oo(edx)
	end % for edx (each frequ ency component)
    end % for cdx (each junction)
end % function coupling_condition

function [AA,bb,ll,cc] = bvp2c_assemble_AA(ypm,odefun,bcfun,xi,xc,dx,neq,nxc,oo,ni,nci,npi)


	% construct function values at section mid-points from separated parts
	yc = zeros(nci(end)-1,1);
	for cdx=1:neq
		switch (oo(cdx))
		case {-1}
			yc(nci(cdx)) =ypm(npi(cdx));
		case {1}
			% first order ode
			%c = cc(1,end-2:end,cdx);
			% TODO this condition has to be checked for each row (!)
			if (false) %0 ~= c(1,2)) 
			yc(nci(cdx):nci(cdx+1)-1) = ( ...
			   		  ypm(npi(cdx)  :2:npi(cdx+1)-2) ...
			   		+ ypm(npi(cdx)+1:2:npi(cdx+1)-1) ...
					);
			else
			if (~dischargeisvariable)
				yc(nci(cdx):nci(cdx+1)-1) = ... 
				   		ypm(npi(cdx)+1:2:npi(cdx+1)-1);
			else
				yc(nci(cdx):nci(cdx+1)-2) = ... 
				   		ypm(npi(cdx)+1:2:npi(cdx+1)-2);
				yc(nci(cdx+1)-1) = ypm(npi(cdx+1)-1);
			end	
			end

		case {2}
			% 2nd order ode
			% sum of left and right going wave as well as inhomogeneous part
			yc(nci(cdx):nci(cdx+1)-1) = ( ...
					  ypm(npi(cdx)  :3:npi(cdx+1)-3) ...
					+ ypm(npi(cdx)+1:3:npi(cdx+1)-2) ...
				 	+ ypm(npi(cdx)+2:3:npi(cdx+1)-1) ...
				 	);
		otherwise
			error('here');
		end	
	end % for cdx

	% ode-coefficients at section mid-points
	cc = feval(odefun,xc,yc);

	% TODO, global buffer
	AA = sparse(npi(end)-1,npi(end)-1,6*npi(end));
	bb = zeros(npi(end)-1,1);
	ll = zeros(nxc,2,neq);

	% for each of the coupled odes
	for cdx=1:neq
		% eigenvalues, roots of the characteristic polynomial
		% roots are in general not complex conjugate pairs,
		% as this is only the case when coefficients are real
		% and case the solution a damped wave
		l            = roots2(cc(:,1:3,cdx));
		ll(:,:,cdx) = l;

		% asumption : ode are either purely 1st or 2nd order, but not mixed
		switch (oo(cdx))
		case {1}
			if (dischargeisvariable)
				[A, b] = bvp1c_assemble_Q(cc,ll,cdx,dx,xi,bcfun,bcarg);
			else
				[A, b] = bvp1c_assemble(cc,ll,cdx,dx,xi,bcfun,bcarg);
			end
		case {2}
			[A, b] = bvp2c_assemble_A(cc,ll,cdx,dx,xi,bcfun,bcarg);
		otherwise
			error('here');
		end
	
		% stack system of odes
		% when the equations are only weekly non-linear,
		% so they can be solved individually
		if (oo(cdx) == -1)
			% TODO deprecated
			AA(npi(cdx), npi(cdx-1):npi(cdx)) = A;
			bb(npi(cdx)) = b;
		else
			AA(npi(cdx):npi(cdx+1)-1, npi(cdx):npi(cdx+1)-1) = A;
			bb(npi(cdx):npi(cdx+1)-1) = b;
		end
	end % for cdx

function [A,b,ih] = bvp2c_assemble_A(cc,ll,ccdx,dx,xi,bcfun,bcarg)
	odec = cc(:,:,ccdx);
	l    = ll(:,:,ccdx);
	m    = 3;
%	odec = bsxfun(@times,odec,1./odec(:,1));
	nobuff = false;
	
	b = zeros(3*nxc,1);
	Abuff = zeros(6+(6+4)*(nxc-1)+3*nxc,3);
	if (nobuff)
	A = sparse([],[],[],m*nxc,m*nxc,(12+2)*nxc);
	end
	% boundary condition at left end
	% f0'(0) = i o z1
	% alternatively, value and derivative could be specified
	% in separate rows, however, than this becomes an ivp and
	% no bc on the right side has to be specified
%	nout = nargout(bcfun);
	% TODO check sign of imaginary part to decide which root is the left going
	nout = 3;
	if (nout < 3)
		[v, p] = bcfun(1,ccdx,bcarg{:});
		%[v, p] = bcfun(xi(1),[],ccdx);
		q = [1, 1];
	else
		[v, p, q] = bcfun(1,ccdx,bcarg{:});
		%[v, p, q] = bcfun(xi(1),[],ccdx);
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end

	if (nobuff)
	A(1,1:m) = [ q(1)*(p(1) + p(2)*l(1,1))*myexp(-0.5*l(1,1)*dx(1)), ...
		     p(1), ...
		     q(2)*(p(1) + p(2)*l(1,2))*myexp(-0.5*l(1,2)*dx(1)), ...
		   ];
	end
	
	Abuf(1:3,1) = 1;
	Abuf(1:3,2) = 1:3;
	Abuf(1,3) = q(1)*(p(1) + p(2)*l(1,1))*myexp(-0.5*l(1,1)*dx(1));
	Abuf(2,3) = p(1);
	Abuf(3,3) = q(2)*(p(1) + p(2)*l(1,2))*myexp(-0.5*l(1,2)*dx(1));
	nbuf = 3;

	% right hand size
	b(1) = v;

	else
		nbuf = 3;
		% dummy zeros
		Abuf(1:3,1)=1;
		Abuf(1:3,2)=1;
	end

	
	k = 1:nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+1;
	Abuf(nbuf+(1:nxc-1),3) = myexp(0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc-1),3) = 1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = myexp(0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -myexp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+5;
	Abuf(nbuf+(1:nxc-1),3) = -1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -myexp(-0.5*l(2:end,2).*dx(2:end));
	nbuf = nbuf+nxc-1;
	
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+1;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,1).*myexp( 0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,2).*myexp( 0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,1).*myexp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,2).*myexp(-0.5*l(2:end,2).*dx(2:end));
	nbuf = nbuf+nxc-1;

	% ode in interior sections
	% at each interior section end point
	%aa = zeros(3,6);
	if (nobuff)
	for k=1:nxc-1
		% continuity of value
		% f_k(k dl) = f_k+1(-k dl)
		% rhs term must be in the centre, as first an last row are overwritten
		% 3,6,9
		A(m*(k-1)+3, m*(k-1) + (1:2*m)) = [	myexp( l(k,1)*dx(k)/2), ...
						        1, ...
						        myexp( l(k,2)*dx(k)/2), ...
			                               -myexp(-l(k+1,1)*dx(k+1)/2), ...
					               -1, ...
					               -myexp(-l(k+1,2)*dx(k+1)/2) ...
						  ];
		%aa(m*(k-1)+3,(1:2*m)) =  A(m*(k-1)+3, m*(k-1) + (1:2*m));
		%if (rescale)
		%	scale = 1./max(abs(A(m*(k-1)+3, m*(k-1) + (1:2*m))));
		%	A(m*(k-1)+3, m*(k-1) + (1:2*m)) = scale*A(m*(k-1)+3, m*(k-1) + (1:2*m));
		%end
		% continuity of first derivative
		% f'_k(k dl) = f'k+1(-k dl)
		%if (rescale)
		%	scale = 1./abs(l(k,1));
		%end
		% 4,7,10
		A(  m*(k-1)+4, m*(k-1) + (1:2*m)) = 1*[ ...
					   l(k,1)*myexp( l(k,1)*dx(k)/2), ...
				         0, ...
					   l(k,2)*myexp( l(k,2)*dx(k)/2), ...
					-l(k+1,1)*myexp(-l(k+1,1)*dx(k+1)/2), ...
					 0, ...
					-l(k+1,2)*myexp(-l(k+1,2)*dx(k+1)/2) ...
					 ];
		%if (rescale)
		%	scale = 1./max(A(  m*(k-1)+4, m*(k-1) + (1:2*m)));
		%	A(  m*(k-1)+4, m*(k-1) + (1:2*m)) = scale.*A(  m*(k-1)+4, m*(k-1) + (1:2*m));
		%end
	end % for k
	end

	% ode at each section mid-point
	% TODO this is only true when c0 != 0
	% general solution : al*myexp(+l*x) + ar*(myexp(+r*x)) + a0 = 0
	%scale_ = [];
	if (nobuff)
	for k=1:nxc
		% 2,5,8
	%	if (rescale)
	%		scale = 1;%./abs(l(k,1).^2.*odec(k,1));
	%	end
		A(m*(k-1)+2, m*(k-1) + 1) = (odec(k,1)*l(k,1).^2 + odec(k,2)*l(k,1) + odec(k,3));
		A(m*(k-1)+2, m*(k-1) + 2) = (odec(k,3));
		A(m*(k-1)+2, m*(k-1) + 3) = (odec(k,1)*l(k,2).^2 + odec(k,2)*l(k,2) + odec(k,3));
		%if (rescale)
			%scale_(k,1) = 1./max(abs(A(m*(k-1)+2, m*(k-1) + (1:3))));
			%A(m*(k-1)+2, m*(k-1) + (1:3)) = A(m*(k-1)+2, m*(k-1) + (1:3)); 
		%end
	end
	end
	k = (1:nxc);
	Abuf(nbuf+(1:nxc),1) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = m*(k-1)+1;
	Abuf(nbuf+(1:nxc),3) = (odec(k,1).*l(k,1).^2 + odec(k,2).*l(k,1) + odec(k,3));
	nbuf = nbuf+nxc;
	Abuf(nbuf+(1:nxc),1) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc),3) = odec(k,3);
	nbuf = nbuf+nxc;
	Abuf(nbuf+(1:nxc),1) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc),3) = (odec(k,1).*l(k,2).^2 + odec(k,2).*l(k,2) + odec(k,3));
	nbuf = nbuf+nxc;

	% rhs
	if (m>2)
		% 2,5,8,...
		%if (rescale)
		%	scale = 1./abs(l(k,1).^2.*odec(k,1));
		%end
		b(m-1:m:m*nxc-1) = -odec(:,4);
	end

	% boundary condition at right end
	% fn(L)  = 0 (or better asymptotic y' = r y)
	if (nout < 3)
		%[v, p] = bcfun(xi(2),[],ccdx);
		[v, p] = bcfun(2,ccdx,bcarg{:});
		q      = [1, 1];
	else
		[v, p, q] = bcfun(2,ccdx,bcarg{:});
		%[v, p, q] = bcfun(xi(2),[],ccdx);
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end
	if (nobuff)
	A(m*nxc,m*nxc+(-m+1:0)) = [ ...
			q(1)*(p(1) + p(2)*l(nxc,1))*myexp(+0.5*l(nxc,1)*dx(nxc)), ...
			p(1), ...
		        q(2)*(p(1) + p(2)*l(nxc,2))*myexp(+0.5*l(nxc,2)*dx(nxc)), ...
			      ];
	end
	%A(m*nxc,m*nxc+(-m:1:0)) = scale*A(m*nxc,m*nxc+(-m:1:0));
	Abuf(nbuf+1,1) = m*nxc;
	Abuf(nbuf+1,2) = m*nxc-m+1;
	Abuf(nbuf+1,3) = q(1)*(p(1) + p(2)*l(nxc,1))*myexp(+0.5*l(nxc,1)*dx(nxc));
	Abuf(nbuf+2,1) = m*nxc;
	Abuf(nbuf+2,2) = m*nxc-m+2;
	Abuf(nbuf+2,3) = p(1); 
	Abuf(nbuf+3,1) = m*nxc;
	Abuf(nbuf+3,2) = m*nxc-m+3;
	Abuf(nbuf+3,3) = q(2)*(p(1) + p(2)*l(nxc,2))*myexp(+0.5*l(nxc,2)*dx(nxc));
	% rhs
	b(m*nxc) = v;
	else
		% dummy zeros
		Abuf(nbuf+(1:3),1) = 1;
		Abuf(nbuf+(1:3),2) = 1;
	end	

	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),3*nxc,3*nxc);
   end % bvp2c_assemble_A
end % bvp2c_assemble_AA

end % bvp2c_assemble

