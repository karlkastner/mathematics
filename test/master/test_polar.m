%perimeter
%timbre

function test_polar()
 opengl neverselect
 N = 2.^(3:8);
 k = 10;
 L0 = 2;

% TODO : not yet second order accurate
% test laplacian only
% V = 

 for ndx=1:length(N)
 	n = N(ndx);
	[L theta r L_] = laplacian_polar(n,L0);
	%A = L;
	V_pot = kron(diag(sparse(1./r)),speye(n)); R_ = []; A_ = 0.5*L_ + V_pot;
	%[V E_] = eigs(A_,R_,k,'SM');
	% coordinate transformed
	V_pot = speye(n^2); R = kron(diag(sparse(r)),speye(n)); A = 0.5*L + V_pot;
	[V E_] = eigs(A,R,k,'SM');

%	P = kron(diag(1./sparse(r)),speye(n))*A;
%	full(P(1:10,1:10))
%	full(A_(1:10,1:10))
%	pause

	[E(:,ndx) ip] = sort(diag(E_));
	V = V(:,ip);

	 X = kron(r,sin(theta)');
	 Y = kron(r,cos(theta)');
	 size(X)
	 figure(ndx);
	 for idx=1:k
		subplot(ceil(sqrt(k)),ceil(sqrt(k)),idx);
		W = reshape(V(:,idx),n,n);
		%imagesc(W);
		contourf(X,Y,W);
		title(num2str(E(idx,ndx)))
		%contourf(X,Y,reshape(diag(R),n,n)); colorbar()
		daspect([1 1 1])
		axis tight
	 end
	E
 end
 Err = E - E(:,end)*ones(1,size(E,2))
 Err = sqrt(sum(Err.*Err))
 figure(ndx+1)
 loglog(N,Err)
 grid on
end % test_polar

