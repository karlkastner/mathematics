% Sat 28 Oct 14:43:06 CEST 2017
%
% function [AA,bb,ll] = bvp2c_assemble(ypm,odefun,bcfun,xi,xc,dx,neq,nxc,oo,ni,npi)
%
%
% output :
%	AA : left hand side discretization matrix
%	bb : right hand side vector
%	ll : eingevalues of the characteristic functions
%
% TODO, treatment of degenerated second order equations a y'' + by' + cy = 0, with c = 0
function [AAA,bbb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni,nci,npi,npii)

	nc = length(odefun);

	% allocate memory for differential operator
	AAA  = sparse([],[],[],npii(end,end)-1,npii(end,end)-1,6*npii(end,end)-1);

	% allocate memory for inhomogeneous part
	bbb  = zeros(npii(end,end)-1,1);

	% for each edge in graph (channel in network)
	for cdx_=1:nc
		% set up discretisation matrix for coupled odes along edge (channel)
		[AA,bb,ll] = bvp2c_assemble_AA(      ypm(npii(1,cdx_):npii(end,cdx_)-1) ...
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
		% TODO global buffers
		AAA(npii(1,cdx_):npii(end,cdx_)-1,npii(1,cdx_):npii(end,cdx_)-1) = AA;
		bbb(npii(1,cdx_):npii(end,cdx_)-1,1)     = bb;
		out(cdx_).ll = ll;
		%.lll(npii(1,cdx_):npii(end,cdx_)-1,1:2,:) = ll;
	end % for cdx_

function apply_coupling_condition()
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
%-> coupling for z0
%-> coupling bc of z0_i == z0_j (n-1 conditions)
%		  sum Q = 0	(nth-condition)
%-> equations for unknown Q0 (nxc+1) :
%-> replace this row by inflow, where Q0 is given
%-> watch out in iteration, initial value of Q0 has to be non-zero



			% mean discharge for each channel, TODO, this has to go into bvp2c_assemble
			for cdx=1:nc
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
		b(rid(1)+(1:3))   = 0; % no external inflow/outflow at bi
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
		end % rdx
		end % if
	end % for idx
    end % for cdx
end % function coupling_condition

function [AA,bb,ll] = bvp2c_assemble_AA(ypm,odefun,bcfun,xi,xc,dx,neq,nxc,oo,ni,nci,npi)

	% ode-coefficients at section mid-points
%	cc = feval(odefun);

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
			if (true) %0 ~= c(1,2)) 
			yc(nci(cdx):nci(cdx+1)-1) = ( ...
			   		  ypm(npi(cdx)  :2:npi(cdx+1)-2) ...
			   		+ ypm(npi(cdx)+1:2:npi(cdx+1)-1) ...
					);
			else
			yc(nci(cdx):nci(cdx+1)-1) = ... 
			   		ypm(npi(cdx)+1:2:npi(cdx+1)-1);
				
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
%	global saveflag
%	if (saveflag)
%		l = load('test.mat');
%		save('test.mat','cc','yc');
%$AA','rr');
%		error('here');
%	end


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
		case {-1}
			[A, b] = bvp_row(cc,ll,cdx,dx,xi,bcfun);
		case {1}
			[A, b] = bvp1c_assemble(cc,ll,cdx,dx,xi,bcfun);
		case {2}
			[A, b] = bvp2c_assemble_A(cc,ll,cdx,dx,xi,bcfun);
		otherwise
			error('here');
		end
	
		% stack system of odes
		% when the equations are only weekly non-linear,
		% so they can be solved individually
		if (oo(cdx) == -1)
			AA(npi(cdx), npi(cdx-1):npi(cdx)) = A;
			bb(npi(cdx)) = b;
		else
			AA(npi(cdx):npi(cdx+1)-1, npi(cdx):npi(cdx+1)-1) = A;
			bb(npi(cdx):npi(cdx+1)-1) = b;
		end
	end % for cdx

function [A,b,ih] = bvp2c_assemble_A(cc,ll,ccdx,dx,xi,bcfun)
	odec = cc(:,:,ccdx);
	l    = ll(:,:,ccdx);
	m    = 3;
%	odec = bsxfun(@times,odec,1./odec(:,1));
	nobuff = false;
	
	b = zeros(nxc,1);
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
		[v, p] = bcfun(xi(1),[],ccdx);
		q = [1, 1];
	else
		[v, p, q] = bcfun(xi(1),[],ccdx);
	end
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
	% b    = -c(1:2:end-1,4);, nope, as this is not the ode, but y1 = y2
	b(1) = v;

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
		[v, p] = bcfun(xi(2),[],ccdx);
		q      = [1, 1];
	else
		[v, p, q] = bcfun(xi(2),[],ccdx);
	end
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
	
	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3));

	% rhs
	b(m*nxc) = v;
	ih = 0;
   end % bvp2c_assemble_A
end % bvp2c_assemble_AA

end % bvp2c_assemble

