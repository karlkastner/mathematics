
N = [ 2.^(1:12)];
k=N(1)-1;
for idx=1:length(N)
	n=N(idx);
	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	tic
	e = eigs(A,k,'SM');
	T(idx,1) = toc;
	Err(idx,1) = (e(1) + pi^2)/pi^2;
	e = eigs(A^2,k,'SM');
	T(idx,2) = toc;
	Err(idx,2) = (e(end) - pi^4)/pi^4;
	e = eigs(A^3,k,'SM');
	T(idx,3) = toc;
	Err(idx,3) = (e(1) + pi^6)/pi^6;
	e = eigs(A^4,k,'SM');
	T(idx,4) = toc;
	Err(idx,4) = (e(end) - pi^8)/pi^8;
	e = eigs(A - 1/12*1/(n+1)^2*A^2,k,'SM')
	T(idx,5) = toc;
	Err(idx,5) = (e(1) + pi^2)/pi^2;
	e = eigs(A - 1/12*1/(n+1)^2*A^2 + 1/90*1/(n+1)^4*A^3,k,'SM')
	T(idx,6) = toc;
	Err(idx,6) = (e(1) + pi^2)/pi^2;
	e = eigs(A - 1/12*1/(n+1)^2*A^2 + 1/90*1/(n+1)^4*A^3 - 1/560*1/(n+1)^6*A^4,k,'SM')
	T(idx,7) = toc;
	Err(idx,7) = (e(1) + pi^2)/pi^2;
end
subplot(2,2,1)
loglog(N,abs(Err));
legend('D2','D2^2','D2^4','D2^6','D4')
ylim([1e-8 1])
subplot(2,2,2)
loglog(N,T);
T./(T(:,1)*ones(1,7))
Err = Err./(Err(:,1)*ones(1,7));
Err(2:4,1:7)
legend('D2','D2^2','D2^4','D2^6','D4')

