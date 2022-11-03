% 2015-08-07 10:25:23.125806143 +0200
% Karl Kastnerm, Berlin

function test_ar1_var_factor()

rho = 0.999;
n = 10;
v = [];
for m=1:n
	v(m,1) = ar1_var_factor(rho,n,m);
	v(m,2) = ar1_var_factor_(rho,n,m);
end
v

nr = 10;
R = (-nr:nr)/(nr+1);
n=10;
v = [];
for idx=1:length(R)
 for m=1:n
	v(m,idx) = ar1_var_factor(R(idx),n,m);
 end
end
plot(v,'.-')
v
end


