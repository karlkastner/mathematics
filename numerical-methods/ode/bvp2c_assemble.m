% Sat 28 Oct 14:43:06 CEST 2017
function [AA,rr,ll] = bvp2c_assemble(ypm,odefun,bcfun,xi,xc,dx,neq,nxc,mm,mi)
myexp=@(x) 1 + x;
	% ode-coefficients at section mid-points
	cc = feval(odefun,xc(1),0);

	% function values at section mid-points
	yc = zeros(neq*nxc,1);
	for idx=1:neq
		% yc = ypm(1:m:end-m+1) + ypm(3:m:end);
		if (2 == mm(idx))
			% first order ode
			c = cc(1,end-2:end,idx);
			if (0 ~= c(1,2)) 
			% not degenerated
			yc((idx-1)*nxc+1:idx*nxc) = ( ...
				  ypm(mi(idx)*nxc+1:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+2:mm(idx):mi(idx+1)*nxc) );
			else
			yc((idx-1)*nxc+1:idx*nxc) = ... 
			   		ypm(mi(idx)*nxc+2:2:mi(idx+1)*nxc);
				
			end
		else
			% 2nd order ode
			% sum of left and right going wave as well as inhomogeneous part
			yc((idx-1)*nxc+1:idx*nxc) = ( ...
				  ypm(mi(idx)*nxc+1:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+2:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+3:mm(idx):mi(idx+1)*nxc) );
		end
	end

%	if (m > 2)
		% inhomogeneous part
		%yi = yc(m:m:end);
%	end
	
	% ode-coefficients at section mid-points
	cc = feval(odefun,xc,yc);

	% TODO allocate AA
	AA = [];
	rr = zeros(mi(end)*nxc,1);
	%ii = []; %zeros(neq*m*nxc,1);
	ll = zeros(nxc,2,neq);

% TODO, global buffer

	% for each of the coupled odes
	for ccdx=1:neq
		% eigenvalues, roots of the characteristic polynomial
		% roots are in general not complex conjugate pairs,
		% as this is only the case when coefficients are real
		% and case the solution a damped wave
		l            = roots2(cc(:,1:3,ccdx));
		ll(:,:,ccdx) = l;

		% silent asumption : ode are either purely 1st or 2nd order, but
		% not mixed
		if (0 == cc(1,1,ccdx))
			[A, b] = bvp1c_assemble(cc,ll,ccdx,dx,xi,bcfun);
		else
			[A, b] = bvp2c_assemble_(cc,ll,ccdx,dx,xi,bcfun);
		end
	
		% stack system of odes
		% when the equations are only weekly non-linear,
		% so they can be solved individually
		AA(1+mi(ccdx)*nxc:mi(ccdx+1)*nxc, ...
		   1+mi(ccdx)*nxc:mi(ccdx+1)*nxc) = A;
		rr(1+mi(ccdx)*nxc:mi(ccdx+1)*nxc,1) = b;
	end % for ccdx

function [A,b,ih] = bvp2c_assemble_(cc,ll,ccdx,dx,xi,bcfun)
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
end % bvp2c_assemble

end % bvp2c_assemble_

