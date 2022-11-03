%-> preselect high frequent part (maxima)
%-> y_hf = -min((dx^2/(a*dt)*I + [1,-2,1])*y,0)
%-> evolve y_hf with euler bw
%-> evolve y_lf with trap
%	-> test and proof that schmeme is positivity preserving
%-> n=10; L = n; Df = fourier_derivative_matrix_1d(n,L,2); x = randn(n,1).^2; dx=ifft(Df*fft(x)); min(dx)
%	-> least squares derivative with wide stencil when not smooth

% why positivity: physically plausible values, stability of relaxation in the reaction equation
% does it also work for advection?

% TODO experiment : accuracy for decreasing dx and dt
%			dirac as initial condition
%			dirac as forcing (point source)
%			dirac as coefficient
%			noise as forcing
% test 2D

% assume a*dt = 1
a=1; dt=1;
n = 1e3;
y = randn(n,1).^2;
I = speye(n);
D2 = derivative_matrix_2_1d(n,n-1,2,'circular');
y_hf = -min(0,(I + 0.5*a*dt*D2)*y);
y_lf = y - y_hf;

y_ = (I - 0.5*a*dt*D2) \ ((I + 0.5*a*dt*D2)*y_lf) + (I - a*dt*D2) \ y_hf;
min(y_)
% more accurate:
y_ = (I - 0.5*a*dt*D2) \ (((I + 0.5*a*dt*D2)*y_lf) + (I - 0.5*a*dt*D2) \ y_hf);
min(y_)
pause


if (1)
figure(1)
n=1e5;
 y = randn(n+1,1).^2;
 x=(0:n)'/n;
K = [0,5];
m = [];
% yi = [interp1(x,y,xi,'spline'),interp1(x,y,xi,'pchip'),interp1(x,y,xi,'linear'),interp1(x,y,xi,'nearest')];
for k=1:20
	D2i = derivative_matrix_2_1d(k*n+1,k*n);
	xi=(0:1/k:n)'/n;
	%yi = [interp1(x,y,xi,'spline'),interp1(x,y,xi,'pchip'),interp1(x,y,xi,'linear'),interp1(x,y,xi,'nearest')];
	yi = [interp1(x,y,xi,'linear')];
	m(k,1) = min(D2i*yi);
		
end
% D2i_ = derivative_matrix_2_1d(k*n+1,k*n,4);
% min(D2i*yi), min(D2i_*yi)
k = 1;
D2i = derivative_matrix_2_1d(k*n+1,k*n);
yi = y;
yi(1) = NaN;
yi(end) = NaN;
yi(:,2) = yi;
yi(:,3) = yi(:,1);
yi(:,4) = yi(:,1);
yi(:,5) = yi(:,1);
w = triwin(1:3);
for idx=1:40;
	m(idx,2:6) = min(D2i*yi)
	yi(:,1) = cvec(conv(yi(:,1),ones(3,1)/3,'same')');
	yi(:,2) = cvec(conv(yi(:,2),w,'same')');
	nf = 2*idx+1;
	yi(:,3) = cvec(conv(y,ones(nf,1)/nf,'same')');
	yi(:,4) = cvec(conv(y,triwin(1:nf),'same')');
	yi(:,5) = cvec(conv(y,hanwin(1:nf),'same')');
end
semilogy(m/m(1),'.')
else
figure(2)
clf
 n=1e2;
 y = randn(n,1).^2;
 dt = 10;
 a=1;
 I = speye(n);
 D2=derivative_matrix_2_1d(n,n+1,2,'circular');
 y_ = (I - 0.5*a*dt*D2) \ ((I + 0.5*a*dt*D2)*y);
res = min(y_,0);
yneg = (I + 0.5*a*dt*D2) \ ((I - 0.5*a*dt*D2)*res);
% res_=(I + 0.5*a*dt*D2) \ ((I - 0.5*a*dt*D2)*yneg);
plot([y,-yneg])
rms([y,yneg])
end
