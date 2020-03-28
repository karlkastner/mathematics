l=10;
k=8;
E = zeros(l,1);
N = zeros(l,1);
R = zeros(l,1);
for idx=1:l
	n=2^idx;
	err = 0;
	r = 0;
	err_l = 0;
%{
	for jdx=1:k
		%A = spdiags(rand(n,3), -1:1, n, n);
		%A = full(A);
		A = rand(n);
		b = rand(n,1);
		%A = 0.5*(A+A');
		%A = sqrt(A'*A);
		
		Ab = A*b;
		err_ = norm(Ab - single(A)*single(b));
		err = err + err_;
		r = r + err_/norm(Ab);
		err_l = err_l + norm(Ab - mmul_accurately(single(A), single(b)));
	end
%}
	A = rand(n);
	b = rand(n,1);
	Ab = b;
	Ab = Ab / norm(Ab);
	Abs = single(b);
	Abs = Abs/norm(Abs);
	As = single(A);
	As = As;
	for jdx=1:k
		Ab = A*Ab;
		Ab = Ab / norm(Ab);
		%Abs = As*Abs;
		Abs = mmul_accurately(As,Abs);
		Abs = Abs / norm(Abs);
	end
	err = norm(Ab - Abs);
	
	N(idx) = n;
	E(idx) = err; %err/k;
	%R(idx) = r/k;
	%R(idx) = err_l/k;
end

clf
loglog(N,[E R sqrt(N)*min(E)]);
loglog(N,[E R log(N)*min(E)]);
%loglog(N,[E R N.^2*min(E)]);
hold on
%loglog(N,min(E)*N.^2,'r');
%loglog(N,sqrt(N)*min(E),'r');
%loglog(N,log(N).^2*min(E),'g');

