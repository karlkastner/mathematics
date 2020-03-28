% 2012-05-29 06:25:55 UTC
% Karl Kästner

% derive finite difference coefficients for a variably spaced grid

clear

syms hk hs hl hr xk xl xc xr xs hrs hkl
x = [xk xl xc xr xs]
x = 1:5;
x =[];

T4 = derive_taylor(4,x)
%T41 = T4(1,:); T42 = T4(2,:); 43 = T4(3,:);T44 = T4(4,:);
T4(2,:) = T4(2,:)/2
T4(3,:) = T4(3,:)/6
for jdx=1:4
	for idx=1:5
		[n d] = numden(T4(jdx,idx));
		N(jdx,idx) = factor(n)
		D(jdx,idx) = factor(d)
	end
end
N(2,:)=N(2,:)*2
N(3,:)=N(3,:)*6
%A = ([                  -hl*hr*hrs,                   -hkl*hr*hrs,         hkl*hrs*(hr - hl) + hl*hr*(hrs-hkl),                   -hkl*hl*hrs,                  -hkl*hl*hr;
% 2*(hl*(hr + hrs) - hr*hrs), 2*(hkl*(hr + hrs) - hr*hrs), 2*(hkl*(hl - hr - hrs) - hl*hr + hrs*(hr -hl)), 2*(hkl*(hl - hrs) - hl*hrs), 2*(hkl*(hl - hr) - hl*hr);
%          -6*(hl - hr - hrs),           -6*(hkl - hr - hrs),                                 6*(hkl + hl - hr - hrs),            6*(hkl + hl - hrs),           6*(hkl + hl - hr);
%                      -24,                        -24,                                                   24,                         24,                       24 ])

pause

P4 = derive_poly(x,x(3))
P40 = P4(1,:);
P41 = P4(2,:);
P42 = P4(3,:);
P43 = P4(4,:);
P44 = P4(5,:);
pause

K4 = derive_power(x)
K41 = K4(1,:);
K42 = K4(1,:);
K43 = K4(1,:);
K44 = K4(1,:);
pause

%K42_ = K2 - 1/12*h^2*A4(3,:)	% only true for constant step-width, and

%K44  = subs(subs(subs(subs(K44 ,hk,xl-xk),hl,xc-xl),hr,xr-xc),hs,xs-xr);
%K44_ = subs(subs(subs(subs(K44_,hk,xl-xk),hl,xc-xl),hr,xr-xc),hs,xs-xr);

%x = 1:5;
%K44  = subs(subs(subs(subs(subs(K44 ,xk,x(1)),xl,x(2)),xc,xc),xr,x(4)),xs,x(5));
%K44_ = subs(subs(subs(subs(subs(K44_,xk,x(1)),xl,x(2)),xc,xc),xr,x(4)),xs,x(5));
%K44.'
%K44_.'
%[num den] =numden(K44_ - K44)
%[num den] = numden((K44_ - K44).')
%pause

%K4.'
%K4_.'

%x = 1:5;
%x = sort(rand(1,5));

%K4  = subs(subs(subs(subs(subs(K4 ,xk,x(1)),xl,x(2)),xc,x(3)),xr,x(4)),xs,x(5))
%K4_ = subs(subs(subs(subs(subs(K4_,xk,x(1)),xl,x(2)),xc,x(3)),xr,x(4)),xs,x(5))
%K4.'
%K4_.'
%K4_ - K4
%[num den] = numden(K4_)
%subs(subs(subs(subs(subs(num,xk,1),xl,2),xc,3),xr,4),xs,5)

