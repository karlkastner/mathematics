% Wed 17 Aug 09:53:26 CEST 2016

- broyden:
	- how to remove the nth-dof at iteration 2n+1?
	- initial matrix:
		chose dx to be eps*1, make H a row vector 
		H = eps*cvec(x)*rvec(f)


% result of analysis:
% broyden's method:
% broydens inverse of the hessian is always far from the true hessian
% linear systems:
%	small:
% 	does not converge within the first 2*npar-1 steps (residual increases)
%	and than finds the solution at step 2*npar (sudden convergence)
%	large:
%	convergence sets in on average at step npar, system solved at 2n
% weakly nonlinear systems:
% 	convergence is super-linear and not sudden
%	the stronger the nonlinearity, the more steps are required
% strongly non-linear: norm of Hi diverges and there is no convergence (restart necessary?)
%
% bfgs: slower convergence (no sudden convergence, better for nonlinear problems)
%
% rectangular systems:
%	hessian inverse of broyden explodes after 2*rank(A) steps, reset?
%	bfgs many iterations required before break down (instant convergence) sets in,
%	almost no convergence before this point, hessian also explodes after 2*n steps
%

x = rand(10,1);

fun = @(c) [x.^4, x.^2, x.^0]*c;
c0 = [1 2 3]';
afun = @(c) fun(c) - fun(c0);

opt.verbose = true;
%c = c0 + 1e-1*c0;
%c = [1 1 1]';
c = zeros(3,1);

%ls_generalized_secant(afun,c,opt)
%ls_broyden(afun,c)
sel = 2;

switch (sel)
case {1}

A = [ 22 10 2 3 7
 14 7 10 0 8
 -1 13 -1 -11 3
 -3 -2 13 -2 4
 9 8 1 -2 4
 9 1 -7 5 -1
 2 -6 6 5 1
 4 5 0 -2 2];
%A = rand(size(A,1),1);
%A = [A(:,1),A(:,1).^2];
%A = A(:,1:3);
A = [A;A.^2];
b = [
 -1 1 0
 2 -1 1
 1 10 11
 4 0 4
 0 -6 -6
 -3 6 3
 1 11 12
 0 -5 -5];
b = [b;b];
b = b(:,1);
afun = @(x) A*x - b
c = pinv(A)*b;
norm(afun(c))
pause
c = zeros(size(A,2),1);
[x nf f H nf_ c_] = ls_broyden(afun,c,0,ns,[],broyden);
nf__ = nf_;
c__ = c_;
H__ = H;
case {2}
nf__ = [];
c__ = [];
H__ = [];
%n = 50;
n  = 4;
nn = 20;
% 5  10 160
% 10 20 400
% 5  20 500
% 4  20 580
% 2 20 200
% 2 40 500
% 3 30
m  = 10; % number of repetitions for statistic
ns = 600; %3*n;
for idx=1:m
if (1)
%n = 6;
%A = 2*(rand(n,3)-0.5); A(:,2) = A(:,1).^2; A(:,3) = A(:,1).^3;
%A = eye(n);
%A = diag([10000 100 1]);
% = diag([16 4 1]);
%A = diag(1:10);
%A = diag(4.^(0:2));
%A = diag(1:3);
%A = eye(10);
%A = diag([1:n]) + 0*randn(n);
%x = linspace(-1,1,n)';
%x = (0:n-1)'/n;
%A = [1./sqrt(2)*x.^0 sin(2*pi*x) cos(2*pi*x)];
%A = [1./sqrt(2)*x.^0 sin(2*pi*x)];
%A = rand(nn,n);
A = repmat(eye(n),nn/n);
%A = eye(nn);
%A = repmat(diag(1:n),nn/n);
A = A + 1e-7/n*rand(size(A));
%A = eye(n) + 1/n*rand(n);
%pause
%A'*A
%pause
%A = vander_1d(x,2);
%A = eye(n);
%A = A + 1/n*randn(n);
% A = A*A';
b  = A*randn(size(A,2),1) + 0.1*randn(size(A,1),1);
c0 = A \ b;
c = 1.1*c0;
end
afun = @(c) A*c + 0.0*A*c.^2 - b;
[q r] = qr(A);
c0 = pinv(r)*q'*b;
c0 = A\b;
r0 = afun(c0);
% for overdetermined (true least squares problem, the optimal residual is nonzero)
n0 = norm(r0);
%norm(A*c0-b)
H0 = A';
H0 = (pinv(A'*A)*A');
H0 = H0  + 1*sqrt(mean(H0(:).^2))*randn(size(H0));
%broyden = false;
broyden = true;
H0 = [];
c = 0.*c0;
[x nf f H nf_ c_ x_] = ls_broyden(afun,c,0,ns,H0,broyden);
[x nf f H nf_ c_ x_] = ls_broyden(afun,x_,0,ns,[],broyden);
%[x nf f H nf_ c_ x_] = ls_broyden(afun,x_,0,ns,[],broyden);
%[x nf f H nf_ c_ x_] = ls_broyden(afun,x_,0,ns,[],broyden);
%$sum([afun(x),afun(c0)].^2)
%pause
%nf_
%n0
%[n0 nf]
%[norm(afun(c0)),norm(afun(x))]
%pause
H__(:,:,idx) = H;
%[x nf f H nf_ c_] = ls_broyden(afun,x,0,30);
%nf_
%n0
%c0
%x
%pause
nf__(:,idx) = nf_-n0;
c__(:,idx) = c_;
%[c0 x]
%nf
%H'*H
%inv(A'*A)
%semilogy(nf_)
end % for
semilogy(mean(nf_,2))

end % case(sel)

figure(1)
subplot(2,2,1)
semilogy((quantile(nf__',[0.16 0.5 0.84]))');
title('Residual norm');
legend('p16', 'median', 'p84');

subplot(2,2,2);
q=quantile(nf__',[0.5])';
plot(q(1:end-1)./q(2:end))
hline(1);
title('Relative reduction of residual norm to previous step');

subplot(2,1,2)
semilogy((quantile(c__',[0.16 0.5 0.84]))')
H=squeeze(nanmedian(H__,3))
title('Norm of inverse Hessian estimate');

