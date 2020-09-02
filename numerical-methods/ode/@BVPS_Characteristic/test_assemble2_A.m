	if (nobuff)
		A = sparse([],[],[],m*nxc,m*nxc,(12+2)*nxc);
	end
		
	if (nobuff)
	A(1,1:m) = [ q(1)*(p(1) + p(2)*l(1,1))*obj.exp(-0.5*l(1,1)*dx(1)), ...
		     p(1), ...
		     q(2)*(p(1) + p(2)*l(1,2))*obj.exp(-0.5*l(1,2)*dx(1)), ...
		   ];
	end
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
	if (nobuff)
	A(m*nxc,m*nxc+(-m+1:0)) = [ ...
			q(1)*(p(1) + p(2)*l(nxc,1))*obj.exp(+0.5*l(nxc,1)*dx(nxc)), ...
			p(1), ...
		        q(2)*(p(1) + p(2)*l(nxc,2))*obj.exp(+0.5*l(nxc,2)*dx(nxc)), ...
			      ];
	end

