% Thu  8 Feb 12:43:59 CET 2024
T = 10;
z0 = 1;
 
I = 1;
A = -1;

dt_ = 2.^(-1:-1:-11);

e = [];
for jdx=1:length(dt_)

dt = dt_(jdx);
nt = T/dt+1;
z=zeros(1,nt);
z(1) = z0;
z(2) = z(1)*exp(A*dt);
t = dt*(0:nt-1)';
for idx=3:nt
	z(idx) = step_bdf2(I,A,z(:,[idx-1,idx-2]),dt);
end
e(jdx,1) = z(end) - z0*exp(A*T);
end
%z= [z',z0*exp(A*t)];
%plot(t,z,'.-');

clf
loglog(dt_/dt_(1),abs(e)/abs(e(1)),'.-')


