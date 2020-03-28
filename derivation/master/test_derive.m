% 2012-03-27 05:51:35 UTC
% Karl Kästner, Berlin

function test_derive()
path(path,'../')

m_=10;
n=1e2;

H = 2.^-(0:m_-1);
Err_sum = zeros(m_,6);

for jdx=1:n
n

%x = sort(pi*rand(1,7)');
x = sort(pi*rand()*(1:7)');
%x = rand()*(1 + rand())*pi*(1:7)';
%x = rand()*pi*(1:7)' + pi*rand();
%x = [1/32 1/16 1/8 1/4 1/2 1/1 2/1]'*0.9*pi + rand();
%x = pi*((1:7).^2'/49*rand());
m = x(4); % 4
x = x - m;
%y = %-sin(m);
%y = m.^4; %exp(m);

%e1 = 4*m.^3;
%e2 = 12*m.^2;
%e3 = 24*m;
%e4 = 24;
e2 = -sin(m);
e4 = sin(m);
%e2 = exp(m);
%e4 = exp(m);
old = 0;
for idx=1:length(H)
	h = H(idx);
	x_ = m + x*h;
	%y = x_.^4;
	y=sin(x_);
	%y=exp(x_);

	D4 = derive_poly(x_(2:6));
	Err(idx,1) = [0 D4(5,:) 0]*y - e4;
	Err(idx,4) = [0 D4(3,:) 0]*y - e2;
%	Err(idx,1) = [ D4(5,:) ]*y - e4;
%	Err(idx,4) = [ D4(3,:) ]*y - e2;

%	[D4(5,:)*y-e4;
%	 D4(3,:)*y-e2]
%Err
%pause

	D4 = derive_taylor(x_(2:6));
	Err(idx,2) = [0 D4(4,:) 0]*y - e4;
	Err(idx,5) = [0 D4(2,:) 0]*y - e2;
%	Err(idx,2) = [ D4(4,:) ]*y - e4;
%	Err(idx,5) = [ D4(2,:) ]*y - e2;
%	D4(2,:)
%	y(3)
%	D4(2,:)*y
%	D4(2,:)*y + y(3)
%pause
D4_ = D4;
%A_=spdiags(ones(7,1)*D4(1,:),-2:2,7,7);
%	[D4 A] = derive_power(vpa(x_(2:6),50));
	[D4 A] = derive_power(x_(2:6));
%	Err(idx,3) = [ D4(4,:) ]*y - e4;
%	Err(idx,6) = [ D4(2,:) ]*y - e2;
%	y = vpa(y,50);
	Err(idx,3) = [0 D4(4,:) 0]*y - e4;
	Err(idx,6) = [0 D4(2,:) 0]*y - e2;

%[0 D4(4,:) 0]*y - old
%old = [0 D4(4,:) 0]*y
%double([[0 D4(4,:) 0]*y e4])
%pause
%	[D4(4,:)*y e4;
%[D4(2,:)*y e2]

%[A D] = laplacian_non_uniform([0; x_; 10]);
%[A_ D_] = laplacian_non_uniform([-10; x_; 20]);
%A_ = D_*A_;
%A = D*A;
%full(A^2 - A_^2)
%A = padarray(A,[1 1]);
%A
%A^2
%b=A*y;
%b2=A*b;
%e = A*y+y;
%R = [ y A*y e A*e A*(A*y)-0*A*e];
%R(3:5,:)
%pause

%pause
%b = [y A*y A*A*y]';
%b(:,3)
%pause
%y(2:4)'
%b=A*y; b(2:4)'
%e=y(2:4)'+b(2:4)'
%e = y+b; e(2:4)'
%c=A*e
%y(3)
%D4_(:,:)*y
% c(2:4)'
%c=A*(-y);
%c(2:4)'
%D4(4,:)*y
%pause
end

Err_sum = Err_sum + Err.^2;

end

Err_sum = sqrt(Err_sum)/sqrt(n);
Err_sum = Err_sum./(ones(m_,1)*Err_sum(1,:));

loglog(1./H,Err_sum)
legend('poly 4th','taylor 4th','power 4th','poly 2nd','taylor 2nd','power 2nd')

end % function test_derive

