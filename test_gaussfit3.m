%Tue 24 Dec 12:34:10 +08 2019
n = 10;
%a = rand(n,1);
%b = rand(n,1);
%c = rand(n,1);
%x = rand(n,3)

%a = 2;
%b = 2;
%c = 4;

%x=1:3;
f = diag(a./(c*sqrt(2*pi)))*exp(-(x-b).^2./(2*c.^2));
f_ =diag(a)*normpdf(x,b,c);

[a_,b_,c_] = gaussfit3(x,f);
[a,b,c]
[a_,b_,c_]
norm(f-f_)
norm(a-a_)
norm(b-b_)
norm(c-c_)
