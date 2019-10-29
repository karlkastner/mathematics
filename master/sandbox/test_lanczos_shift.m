% 2013 Mar 11 02:45 Karl Kästner

% convergence of Lanczos on inverserely shifted matrix
% example : smallest eigenvalue of disrcrete Laplacian

function test_lanczos_shift()
	n = 100;
	k = 5*n;

	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	I = speye(n);

%	A = kron(A,I) + kron(I,A);
%	I = kron(I,I);
	q1 = rand(length(A),1);

	% common lanczos
	T1    = lanczos(A,I,q1,k);
	e(1)  = max(eig(full(T1)));

	% inverted lanczos
	T2    = lanczos(-I,-A,q1,k);
	e(2)  = 1./min(eig(full(T2)));

	% inner - outer scheme
	l = 4;
	[e(3) e3] = outer(A,I,q1,k,l);

	
	e_true = eigs(A,[],1,'SM');

	% convergence of unmodified scheme
	for idx=1:length(T1)
		e1(idx,1) = max(eig(T1(1:idx,1:idx)));
	end

	% convergence of inverted scheme
	for idx=1:length(T2)
		e2(idx,1) = max(1./eig(T2(1:idx,1:idx)));
	end
	size(e2)


	err1 = e1 - e_true;
	err2 = e2 - e_true;
	err3 = e3 - e_true;
	err1 = err1 / err1(1);
	err2 = err2 / err2(1);
	err3 = err3 / err3(1);
	err3 = err3(1:k);
	size(err1)
	size(err2)
	size(err3)
	semilogy(abs([err1 err2 err3]));
	[e e_true]
	legend('plain', 'inverted', 'restarted');
	% + pi^2
end

function [e e_] = outer(A,B,q1,k,l)
%	k_ = k/l;
%	for jdx=0:l-1
	m = 0;
	k_ = 1;
	while(m < k)
		mu = q1'*A*q1 / (q1'*B*q1);
	c = 1;
	Di = c*B;
	D = 1/c*B;
%
		[T Q] = lanczos(Di*(A-mu*B),D,q1,k_);
		T = T+mu*speye(length(T));
		[v e] = eigs(T,[],1,'LA');
		q1 = Q*v;
		for idx=1:length(T)
			m=m+1;
			e_(m,1)= max(eig(T(1:idx,1:idx)));
		end
		%k_ = k_+1; % linear
		%k_ = k_+k_; % quadratic
		k_ = 2*k_; % exponential
	end
end


function [T Q] = lanczos(A,B,q1,k)
	tol = 1e-15;
	buf = zeros(0,3);
	q1 = q1/sqrt(q1'*B*q1);
	Bq1 = B*q1; 
	Bq2 = A*q1;
	Q = q1;
	for idx=1:k
		alpha = q1'*Bq2;
		Bq2   = Bq2 - alpha*Bq1;
		q2   = B \ Bq2;
		beta  = sqrt(q2'*Bq2);
		buf(end+1,:) = [idx idx alpha];
		if (beta > tol && idx < k)
			buf (end+1,:) = [ idx idx+1 beta];
			buf (end+1,:) = [ idx+1 idx beta];
			q2 = (1/beta)*q2;
			Bq2 = (1/beta)*Bq2;
			Q(:,end+1) = q2;
			Bq0 = Bq1;
			Bq1 = Bq2;
			q0 = q1;
			q1 = q2;
			Bq2   = A*q1 - beta*Bq0;
		else
			break;
		end
%		Q'*B*Q
%		Q'*Q
%		pause
	end % for idx
	T = sparse(buf(:,1),buf(:,2),buf(:,3));
end % function lanczos


function [T Q] = lanczos_(A,B,q1,k)
	tol = 1e-15;
	buf = zeros(0,3);
	q1 = q1/sqrt(q1'*B*q1);
	q2 = A*(B*q1);
	Q = q1;
	for idx=1:k
		alpha = q2'*B*q1;
		q2    = q2 - alpha*q1;
		q2    = B \ q2;
		beta  = sqrt(q2'*B*q2);
		buf(end+1,:) = [idx idx alpha];
		if (beta > tol && idx < k)
			buf (end+1,:) = [ idx idx+1 beta];
			buf (end+1,:) = [ idx+1 idx beta];
			q2 = (1/beta)*q2;
			Q(:,end+1) = q2;
			q0 = q1;
			q1 = q2;
			q2   = A*(B*q1) - beta*q0; % no b here
		end % if
	end % for idx
	T = sparse(buf(:,1),buf(:,2),buf(:,3));
	Q'*B*Q
	pause
end % function lanczos

