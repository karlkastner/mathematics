% the only way to reduce the error of the dimensional splitting is to reduce
% the time step

% maybe this can be improved by combining with the diagonal, 
% but best is to go 2D(!)

n = 50;
r_ = 10;
yy = [];
clf;
for idx=1:2;
	d2 = derivative_matrix_2_1d(n,n,2*idx,'circular');
	I = speye(n);
	x = zeros(n,1);
	x(end/2) = 1;
	A = (I - r_/jdx*d2);
	y = A \ x;
	yy(:,:,idx) = y*y';
%	for jdx=1:2;
%	A = (I - 1/jdx*d2);
%	y(:,jdx) = A^jdx\x;
%	 end;
end

[dx,dy,dxx,dxy,dyy] =derivative_matrix_2d([n,n],[n,n],4,{'circular','circular'});
I = speye(n*n);
x = zeros(n,n);
x(end/2,end/2) = 1;
x=x(:);
q=1;
A=(I-r_/q*dxx-r_/q*dyy)^q;
yy(:,:,3) = reshape((A\x),n,n);
q=4;
A=(I-r_/q*dxx-r_/q*dyy)^q;
yy(:,:,4) = reshape((A\x),n,n);

for idx=1:size(yy,3)
	subplot(3,4,idx);
	imagesc(yy(:,:,idx));
	r=(yy(:,:,idx)-yy(:,:,4));
	[rms(r(:)),max(r(:)),min(r(:))]
	imagesc(r);
	subplot(3,4,idx+4);
	cla
	k = 1:n; k=k-k(end/2);
	k_ = k*sqrt(2);
	plot(k_,diag(yy(:,:,idx))); %/max(flat(yy(:,:,idx))))
	hold on
	plot(k,yy(:,end/2,idx));
	subplot(3,4,idx+8);
	plot(k_,diag(r(:,:,1))); %/max(flat(yy(:,:,idx))))
	hold on
	plot(k,r(:,end/2,1));
end


if (0)
% not more accurate and not pos preserving
p=1;
r=0.2;
 n = 21;
 [d1x,d1y,dxx,dxy,dyy] = derivative_matrix_2d([n,n],[1,1],2,{'circular','circular'});
 I = speye(n^2);
 x=zeros(n,n);
 x(round(end/2),round(end/2)) = 1;

 subplot(2,3,1);
if (0)
 y0 = reshape(     ...
	 	(I + (I - r*dxx - r*dyy)^-1*(r/p)^2*dxx*dyy) ...
		* ...
		   ( (I-(r/p)*dxx) ...
	         \ (   (I-(r/p)*dyy) ...
			\x(:))) ...
                 ,n,n);
	 	%(I + (I - r*dxx)^-1*(I - r*dyy)^-1*(r/p)^2*dxx*dyy) ...
 y1 = reshape(     ...
	 	(I + (r/p)^2*dxx*dyy) ...
		* ...
		   ( (I-(r/p)*dxx) ...
	         \ (   (I-(r/p)*dyy) ...
			\x(:))) ...
                 ,n,n);
end
%y1  = 2*( (I-1/2*r*dxx)^2 \ ((I - 1/2*r*dyy)^2 \ x(:))) - ((I-r*dxx) \ ((I - r*dyy) \ x(:)));
q = 2;
y1 = x(:);
for idx=1:q
	y1  = (I-1/q*r*dxx) \ ((I - 1/q*r*dyy) \ y1(:));
end
%y1  = ((I-1/q*r*dxx)*(I - 1/q*r*dyy))^q \ x(:);
%y1  = ( (I-1/2*r*dxx)^2 \ ((I - 1/2*r*dyy) \ x(:))) - ((I-r*dxx) \ ((I - r*dyy) \ x(:)));
y1 = reshape(y1,n,n);
 imagesc(y1);
colorbar
 axis square;

		 %((I+(r/p)^2*dxy) ...

 subplot(2,3,2);
y2=reshape(   (I-(r/p)*dxx) ...
	    \ ( (I-(r/p)*dyy) ...
	        \ x(:) ...
              ) ...
                 ,n,n);
 imagesc(y2);
colorbar
 axis square;

 subplot(2,3,3);
y3 = reshape( (I-(r/p)*dxx-(r/p)*dyy)^p\x(:),n,n);
q=10;
y3 = x(:);
for idx=1:q
	y3  = (I-1/q*r*dxx) \ ((I - 1/q*r*dyy) \ y3(:));
	%y3 = (I-0.5*(r/q)*dxx-0.5*(r/q)*dyy) \ ((I+0.5*(r/q)*dxx+0.5*(r/q)*dyy)*x(:));
end
y3 = reshape(y3,n,n);
imagesc(y3)
colorbar
axis equal;

%rms(flat(y0-y3))
rms(flat(y1-y2))
rms(flat(y2-y3))
rms(flat(y1-y3))
end

