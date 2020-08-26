% Sat 28 Oct 14:43:06 CEST 2017
% TODO, treatment of degenerated second order equations a y'' + by' + cy = 0, with c = 0

function [A,b] = assemble2_A(obj,cdx,edx)

%	xi     = obj.xi(cdx,:);
%	xc     = obj.out(cdx).xc;
	dx     = obj.out(cdx).dx;
	nxc    = obj.nxc(cdx);

	odec = obj.out(cdx).cc(:,:,edx);
	l    = obj.out(cdx).ll(:,:,edx);
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
% TODO pass channel id instead of indexing
		[v, p] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});
		%[v, p] = bcfun(xi(1),[],ccdx);
		q = [1, 1];
	else
		[v, p, q] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});
		%[v, p, q] = bcfun(xi(1),[],ccdx);
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end

	if (nobuff)
	A(1,1:m) = [ q(1)*(p(1) + p(2)*l(1,1))*obj.exp(-0.5*l(1,1)*dx(1)), ...
		     p(1), ...
		     q(2)*(p(1) + p(2)*l(1,2))*obj.exp(-0.5*l(1,2)*dx(1)), ...
		   ];
	end
	
	Abuf(1:3,1) = 1;
	Abuf(1:3,2) = 1:3;
	Abuf(1,3) = q(1)*(p(1) + p(2)*l(1,1))*obj.exp(-0.5*l(1,1)*dx(1));
	Abuf(2,3) = p(1);
	Abuf(3,3) = q(2)*(p(1) + p(2)*l(1,2))*obj.exp(-0.5*l(1,2)*dx(1));
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
	Abuf(nbuf+(1:nxc-1),3) = obj.exp(0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+2;
	Abuf(nbuf+(1:nxc-1),3) = 1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = obj.exp(0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -obj.exp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+5;
	Abuf(nbuf+(1:nxc-1),3) = -1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -obj.exp(-0.5*l(2:end,2).*dx(2:end));
	nbuf = nbuf+nxc-1;
	
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+1;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,1).*obj.exp( 0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,2).*obj.exp( 0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,1).*obj.exp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = m*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = m*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,2).*obj.exp(-0.5*l(2:end,2).*dx(2:end));
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
		A(m*(k-1)+3, m*(k-1) + (1:2*m)) = [	obj.exp( l(k,1)*dx(k)/2), ...
						        1, ...
						        obj.exp( l(k,2)*dx(k)/2), ...
			                               -obj.exp(-l(k+1,1)*dx(k+1)/2), ...
					               -1, ...
					               -obj.exp(-l(k+1,2)*dx(k+1)/2) ...
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
					   l(k,1)*obj.exp( l(k,1)*dx(k)/2), ...
				         0, ...
					   l(k,2)*obj.exp( l(k,2)*dx(k)/2), ...
					-l(k+1,1)*obj.exp(-l(k+1,1)*dx(k+1)/2), ...
					 0, ...
					-l(k+1,2)*obj.exp(-l(k+1,2)*dx(k+1)/2) ...
					 ];
		%if (rescale)
		%	scale = 1./max(A(  m*(k-1)+4, m*(k-1) + (1:2*m)));
		%	A(  m*(k-1)+4, m*(k-1) + (1:2*m)) = scale.*A(  m*(k-1)+4, m*(k-1) + (1:2*m));
		%end
	end % for k
	end

	% ode at each section mid-point
	% TODO this is only true when c0 != 0
	% general solution : al*obj.exp(+l*x) + ar*(obj.exp(+r*x)) + a0 = 0
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
		[v, p] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});
		q      = [1, 1];
	else
		[v, p, q] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});
		%[v, p, q] = bcfun(xi(2),[],ccdx);
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end
	if (nobuff)
	A(m*nxc,m*nxc+(-m+1:0)) = [ ...
			q(1)*(p(1) + p(2)*l(nxc,1))*obj.exp(+0.5*l(nxc,1)*dx(nxc)), ...
			p(1), ...
		        q(2)*(p(1) + p(2)*l(nxc,2))*obj.exp(+0.5*l(nxc,2)*dx(nxc)), ...
			      ];
	end
	%A(m*nxc,m*nxc+(-m:1:0)) = scale*A(m*nxc,m*nxc+(-m:1:0));
	Abuf(nbuf+1,1) = m*nxc;
	Abuf(nbuf+1,2) = m*nxc-m+1;
	Abuf(nbuf+1,3) = q(1)*(p(1) + p(2)*l(nxc,1))*obj.exp(+0.5*l(nxc,1)*dx(nxc));
	Abuf(nbuf+2,1) = m*nxc;
	Abuf(nbuf+2,2) = m*nxc-m+2;
	Abuf(nbuf+2,3) = p(1); 
	Abuf(nbuf+3,1) = m*nxc;
	Abuf(nbuf+3,2) = m*nxc-m+3;
	Abuf(nbuf+3,3) = q(2)*(p(1) + p(2)*l(nxc,2))*obj.exp(+0.5*l(nxc,2)*dx(nxc));
	% rhs
	b(m*nxc) = v;
	else
		% dummy zeros
		Abuf(nbuf+(1:3),1) = 1;
		Abuf(nbuf+(1:3),2) = 1;
	end	

	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),3*nxc,3*nxc);
end % assemble2_A

