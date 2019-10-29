function [a b T] = lanczos(A)
	n = size(A,1)
	a = zeros(n,1);
	b = zeros(n+1,1);
	q_1 = rand(n,1);
	q_1 = q_1/norm(q_1);
	b(1) = 0;
	q_2 = A*q_1;
	for idx=1:n
		a(idx) = q_2'*q_1;
		q_2 = q_2 - a(idx)*q_1;
		b(idx+1) = norm(q_2);
		q_2 = (1/b(idx+1))*q_2;
		q_0 = q_1;
		q_1 = q_2;
		q_2 = A*q_1 - b(idx+1)*q_0;
	end
	b
	T = spdiags([b(2:end,1) a b(1:end-1,1)],-1:1,n,n);
end

