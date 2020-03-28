
% fehlerentwicklung proportional zur matrixnummer und condition number
% AQ_n = Q_n+1*H_n holds always, but:
% A = Q_n+1*H_n*Q_n does not hold, due to round-off errors

function test_arnoldi(n_max)

n=10;
A0 = rand(n);
A = sqrt(A0'*A0);
[H Q] = arnoldi(A, n+1, 0);
% weak test
norm(A*Q(1:n,1:n) - Q*H)
% strong test
norm(H(1:n,:) - Q(:,1:n)'*A*Q(:,1:n))

%[H Q] = arnoldi(A, n/2, 0);
%norm(A*Q(:,1:n/2-1) - Q*H)

% make spd
%A0= A0+diag(1:length(A0));
A = sqrt(A0'*A0);
[H Q] = arnoldi(A, n+1, 1);
norm(A*Q(1:n,1:n) - Q*H)
norm(H(1:n,:) - Q(:,1:n)'*A*Q(:,1:n))
%[H Q] = arnoldi(A, n+1, 0)
%[H Q] = arnoldi(A, n+1, 1, 0)
%[H Q] = arnoldi(A, n+1, 1, 1)
%pause

% make snd
A = -sqrt(A0'*A0);
[H Q] = arnoldi(A, n+1, 1);
norm(A*Q(1:n,1:n) - Q*H)
norm(H(1:n,:) - Q(:,1:n)'*A*Q(:,1:n))

% make symmetric-mixed-definite
A = sqrt(A0'*A0) - eye(n);
[H Q] = arnoldi(A, n+1, 1);
norm(A*Q(1:n,1:n) - Q*H)
norm(H(1:n,:) - Q(:,1:n)'*A*Q(:,1:n))
	

Err_weak = zeros(n_max,3);
Err_strong = zeros(n_max,3);
Err_ortho = zeros(n_max,3);
% repeat for higher accuracy
R = [];
r_max=1;
for rdx=1:r_max;
for ndx=1:n_max
	n = 2^ndx
	N(ndx) = n;
	% spd-random matrix
	A = rand(n);
	%A = sqrt(A*A');
	A = sqrt(A*A');
	%nor = norm(A);
	%A = diag(1:n) + diag(1:n-1,-1) + diag(1:n-1,+1) + ones(n)

	E = sort(eig(A));
	for idx=0:3
			m = n+1;
		tic;
		switch(idx)
			case {0}
				% arnoldi method
				
				[H Q] = arnoldi(A, m, 0, 0);
				e1 = sort(eig(A));
			case {1}
				% lanczos method
				[H Q] = arnoldi(A, m, 1, 0);
				e2 = sort(eig(A));
			case {2}
				% lanczos with full reorthogonalisation
				[H Q] = arnoldi(A, m, 1, 1);
				e3 = sort(eig(A));
			case {3}
				% lanczos with periodic reorthogonalisation
				[H Q] = arnoldi(A, m, 1, 2);
				e4 = sort(eig(A));
		end
		T(ndx,idx+1) = toc;
		%R(ndx, idx+1) = rank(H(1:n,:))-n;
		Err_weak(ndx,idx+1) = norm(A*Q(1:n,1:n) - Q*H);
		Err_strong(ndx,idx+1) = norm(H(1:n,:) - Q(:,1:n)'*A*Q(:,1:n));
		Err_eig(ndx,idx+1) = norm(sort(eig(H(1:n,:))) - E);
		Err_ortho(ndx,idx+1) = norm(Q(1:n,1:n)'*Q(1:n,1:n) - eye(n));
	end
end
end % rdx

subplot(1,5,1)
loglog(N, [Err_weak]); %, Err_strong, Err_ortho])
legend('arnoldi', 'lanczos', 'lre', 'lpr'); % 'arnoldi strong', 'lanczos strong', 'lre-strong', 'eoa', 'eol', 'eolro')
grid on
set(gca, 'minorgrid','none')
subplot(1,5,2)
loglog(N, [Err_strong]); %, Err_strong, Err_ortho])
grid on
set(gca, 'minorgrid','none')
subplot(1,5,3)
loglog(N, [Err_ortho]); %, Err_strong, Err_ortho])
%semilogx(N,R)
grid on
set(gca, 'minorgrid','none')

subplot(1,5,4)
loglog(N, [Err_eig]); %, Err_strong, Err_ortho])
%semilogx(N,R)
grid on
set(gca, 'minorgrid','none')

subplot(1,5,5)
loglog(N, [T]); %, Err_strong, Err_ortho])
%semilogx(N,R)
grid on
set(gca, 'minorgrid','none')

end % test_arnoldi 

