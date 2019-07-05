% Sat Feb 19 19:59:39 MSK 2011
% Karl Kästner, KTH Stockholm
%%
%% apply boundary conditions
function q = apply_bc(obj,t,q,dt)
	if (SWE.BCVER<1)
		q = feval(obj.bcfun, t, q);
	else
		m   = obj.m;
		n   = length(obj.x);
		dx  = obj.dx;

%	(q(n-1)-q(n-2))/obj.dx
%	q(2*n-1).^2*2e-3/9.81*q(n-1).^-3

		% left end
		% value of first interior point
		ql1 = q(2:n:n*(m-1)+2);
		[A B rhs] = feval(obj.bcfun{1},obj.pde,t,ql1,dt);
		[A B]     = bc_transform(A,B);

		% value of left boundary point
		ql0 = A \ (rhs + B*ql1);
		%ql0 = (rhs - B.*ql1)./A;

		% right end
		qr1 = q(n-1:n:n*m-1);
		[A B rhs] = feval(obj.bcfun{2},obj.pde,t,qr1,dt);
		[A B]  = bc_transform(A,B);
		qr0 = A \ (rhs + B*qr1);
		%qr0 = (rhs - B.*qr1)./A;


		% write back
		q(1:n:n*(m-1)+1) = ql0;
		q(n:n:n*m) = qr0;

%			bl  = [qold(1),   q(2), qold(n+1), q(n+2),   1].';
%			br  = [qold(n), q(n-1), qold(2*n), q(2*n-1), 1].';
%			ql  = L*bl;
%			qr  = R*br;
%			q(1)   = ql(1);
%			q(n+1) = ql(2);
%			q(n)   = qr(1);
%			q(2*n) = qr(2);
%			% right side
%			[p rhs]	= feval(obj.bcfun{2},t,[q(n);q(2*n)]);
%			%fun = obj.bcfun(2);
%			%[p rhs]	= fun(t,[q(n);q(2*n)]);
%			% [-p(2)/dx, p(1)+p(2)/dx] = rhs
%			q(n)   = (rhs(1)-(-p(1,2)/dx)*q(n-1))/(p(1,1)+p(1,2)/dx);
%			q(2*n) = (rhs(2)-(-p(2,2)/dx)*q(2*n-1))/(p(2,1)+p(2,2)/dx);
	end % else of if 0,1

	% transform boundary condistion specified as values and derivatives
	% into factors for of the first two independent variables
	% in:	rhs  = Ain*q|0 + Bin*dq/dx|0
	% out:	Ao*[q0|0; ..; qm|0] = Boq*[q0|1; ..; qm|1] + rhs
	%
%	% in:  Ai  = [q0|_0, .. , q_m|_0]
%	%      Bi  = [dq0/dx|_0, .., qm|_0,dqn/dx|_0] = rhs
%	% out: A [q00; ...; qn0] = rhs - B[q01; ...; qn1];
%	% A and B are diagonal matrices and stored as vectors
%	% in:  Ai  = [q0|_0,dq0/dx|_0; ...; qn|_0,dqn/dx|_0] = rhs
%	% out: A [q00; ...; qn0] = rhs - B[q01; ...; qn1];
%	% A and B are diagonal matrices and stored as vectors
	function [Aout Bout] = bc_transform(Ain,Bin)
		dx   = obj.dx;
		if (1)
		Aout = Ain + 1/dx*Bin;
		Bout = 1/dx*Bin;
%pause
		else
			A = Ain(:,1) - Ain(:,2)/dx;
			B = Ain(:,2)/dx;
		end
	end % bc_transform

end % apply_bc

