% Tue Jun  5 19:00:22 MSK 2012
% Karl Kästner, Berlin

function test_qr()
	path(path,'/home/pia/Desktop/cosse/kth/dn2230-numalg/hw1/');
	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	E = eig(A);
	E_ = eig_symmetric(A,1);
	norm(sor(E) - sort(E_))
end


