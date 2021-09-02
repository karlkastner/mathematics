% Sat 28 Oct 14:43:06 CEST 2017
%
%% assemble the discretisation matrix for a second-order ode
%% (non-zero frequency component)
%
% TODO treatment of degenerated second order equations a y'' + by' + cy = 0, with c = 0
% TODO during iteration, only variable elements of the buffer have to be rewritten
function assemble2_A(obj,cdx,edx)

	dx     = obj.out(cdx).dx;
	nxc    = obj.nxc(cdx);

	odec = obj.out(cdx).cc(:,:,edx);
	l    = obj.out(cdx).ll(:,:,edx);
	
	% TODO directly write to global buffers
	b = zeros(3*nxc,1);
	Abuf = zeros(6+(6+4)*(nxc-1)+3*nxc,3);

	% boundary condition at left end
	% f0'(0) = i o z1
	% alternatively, value and derivative could be specified
	% in separate rows, however, than this becomes an ivp and
	% no bc on the right side has to be specified
%	nout = nargout(bcfun);
	% TODO check sign of imaginary part to decide which root is the left going
	nout = 3;
	if (nout < 3)
		[v, p] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});
		q = [1, 1];
	else
		[v, p, q] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
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


	% continuity of f at segment interfaces	
	k = 1:nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+1;
	Abuf(nbuf+(1:nxc-1),3) = obj.exp(0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+2;
	Abuf(nbuf+(1:nxc-1),3) = 1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = obj.exp(0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -obj.exp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+5;
	Abuf(nbuf+(1:nxc-1),3) = -1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -obj.exp(-0.5*l(2:end,2).*dx(2:end));
	nbuf = nbuf+nxc-1;
	
	% continuity of f' at sedment interfaces
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+1;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,1).*obj.exp( 0.5*l(1:end-1,1).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc-1),3) = l(1:end-1,2).*obj.exp( 0.5*l(1:end-1,2).*dx(1:end-1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,1).*obj.exp(-0.5*l(2:end,1).*dx(2:end));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+(1:nxc-1),1) = 3*(k-1)+4;
	Abuf(nbuf+(1:nxc-1),2) = 3*(k-1)+6;
	Abuf(nbuf+(1:nxc-1),3) = -l(2:end,2).*obj.exp(-0.5*l(2:end,2).*dx(2:end));
	nbuf = nbuf+nxc-1;

	% ode at section mid-points
	% TODO this is only true when c0 != 0
	% general solution : al*obj.exp(+l*x) + ar*(obj.exp(+r*x)) + a0 = 0
	k = (1:nxc);
	Abuf(nbuf+(1:nxc),1) = 3*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = 3*(k-1)+1;
	Abuf(nbuf+(1:nxc),3) = (odec(k,1).*l(k,1).^2 + odec(k,2).*l(k,1) + odec(k,3));
	nbuf = nbuf+nxc;
	Abuf(nbuf+(1:nxc),1) = 3*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = 3*(k-1)+2;
	Abuf(nbuf+(1:nxc),3) = odec(k,3);
	nbuf = nbuf+nxc;
	Abuf(nbuf+(1:nxc),1) = 3*(k-1)+2;
	Abuf(nbuf+(1:nxc),2) = 3*(k-1)+3;
	Abuf(nbuf+(1:nxc),3) = (odec(k,1).*l(k,2).^2 + odec(k,2).*l(k,2) + odec(k,3));
	nbuf = nbuf+nxc;

	% rhs
	if (3>2)
		% 2,5,8,...
		b(2:3:3*nxc-1) = -odec(:,4);
	end

	% boundary condition at right end
	% fn(L)  = 0 (or better asymptotic y' = r y)
	if (nout < 3)
		[v, p] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});
		q      = [1, 1];
	else
		[v, p, q] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});
	end
	if (~isempty(v))
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end
%	A = sparse(Abuf(1:nbuf,1),Abuf(1:nbuf,2),Abuf(1:nbuf,3),3*nxc,3*nxc);
%	A(end,:)
%pause
	Abuf(nbuf+1,1) = 3*nxc;
	Abuf(nbuf+1,2) = 3*nxc-3+1;
	Abuf(nbuf+1,3) = q(1)*(p(1) + p(2)*l(nxc,1))*obj.exp(+0.5*l(nxc,1)*dx(nxc));
	Abuf(nbuf+2,1) = 3*nxc;
	Abuf(nbuf+2,2) = 3*nxc-3+2;
	Abuf(nbuf+2,3) = p(1); 
	Abuf(nbuf+3,1) = 3*nxc;
	Abuf(nbuf+3,2) = 3*nxc-3+3;
	Abuf(nbuf+3,3) = q(2)*(p(1) + p(2)*l(nxc,2))*obj.exp(+0.5*l(nxc,2)*dx(nxc));
	% rhs
	b(3*nxc) = v;
	else
		% dummy zeros
		Abuf(nbuf+(1:3),1) = 1;
		Abuf(nbuf+(1:3),2) = 1;
	end	
	nbuf = nbuf+3;

	Abuf(:,1:2) = Abuf(:,1:2) + obj.npi(edx,cdx)-1 + obj.npii(1,cdx)-1;
	obj.Abuf   = [obj.Abuf; Abuf];
	obj.b(obj.npii(edx,cdx):obj.npii(edx+1,cdx)-1,1)   = b;
%	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),3*nxc,3*nxc);
end % assemble2_A

