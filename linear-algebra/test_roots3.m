% Fri 27 Mar 10:33:15 +08 2020
n = 10;
A = rand(n,4);

r = [];
for idx=1:n
	r(idx,:) = roots(A(idx,:));
end

rr = roots3(A);

r = sort(r.').'
rr = sort(rr.').'
r - rr
%x=(1:4); r = roots(x).', c=roots2poly(r)

AA = roots2poly(rr)

res = AA - A./A(:,1)
rms(res)

