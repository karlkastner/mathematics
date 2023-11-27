% 2023-03-31 10:36:27.729311851 +0200
% this can be reduced by exploiting stationarity
% sum_i sum_j

n = 2;

x1 = 0;
x2 = 0;
dx = 1;
a  = 1;

dx_ = dx*innerspace(-0.5,0.5,n);
cfun = @(x) exp(-a*abs(x));

% so the min difference is j0 (idx==jdx+j0), the max difference is n-1 (idx=i0,jdx=j0+n)
c = 0;
for idx=1:n
	for jdx=1:n
		c = c + cfun((x2+dx_(idx))-(x1+dx_(jdx)));
	end
end
c  = c/n^2;
c_ = c;

% this holds for even and odd
c = n*cfun(abs(x2-x1));
dx__ = dx*(1:n)/n;
k = n;
for idx=1:n-1
	k=k+2*(n-idx);	
	c = c + (n-idx)*cfun(abs(x2-x1)-(dx__(idx)));
	c = c + (n-idx)*cfun(abs(x2-x1)+(dx__(idx)));
end
[k,n^2]
c = c/n^2;
c
c_

n=2;
cfun2 = @(x,y) exp(-a*hypot(x,y));
y1 = 0;
y2 = 0;
dy = dx;
dy_ = dx_;

x1 = (0:3)';
x1 = repmat(x1,1,4);
y1 = (0:3);
y1 = repmat(y1,4,1);

n = 10;
x1 = 0; x2=0;
y1=0; y2=0;

tic()
c_ = cov_cell_averages_2d(cfun2,x2-x1,y2-y1,dx,dy,n,0)
toc()
tic()
c = cov_cell_averages_2d(cfun2,x2-x1,y2-y1,dx,dy,n,1)
toc()

if (0)
% 2d
c = 0;
for idx=1:n
	for jdx=1:n
		for i2=1:n
			for j2=1:n
				c = c + cfun2((x2+dx_(idx))-(x1+dx_(jdx)), ...
					      (y2+dy_(i2))-(y1+dy_(j2)));
			end
		end
	end
end
c  = c/n^4;
c_ = c;



c = n^2*cfun2(x2-x1,y2-y1);
dx__ = dx*(1:n)/n;
dy__ = dy*(1:n)/n;
k = n^2;
for idx=1:n-1
	k = k+4*(n-idx)*n;
	c = c + (n-idx)*n*cfun2((x2-x1)-(dx__(idx)),(y2-y1));
	c = c + (n-idx)*n*cfun2((x2-x1)+(dx__(idx)),(y2-y1));
	c = c + (n-idx)*n*cfun2((x2-x1),(y2-y1)-(dy__(idx)));
	c = c + (n-idx)*n*cfun2((x2-x1),(y2-y1)+(dy__(idx)));
	for jdx=1:n-1
		k = k+4*(n-idx)*(n-jdx);
		c = c + (n-idx)*(n-jdx)*cfun2((x2-x1)-(dx__(idx)),(y2-y1)-(dy__(idx)));
		c = c + (n-idx)*(n-jdx)*cfun2((x2-x1)-(dx__(idx)),(y2-y1)+(dy__(idx)));
		c = c + (n-idx)*(n-jdx)*cfun2((x2-x1)+(dx__(idx)),(y2-y1)-(dy__(idx)));
		c = c + (n-idx)*(n-jdx)*cfun2((x2-x1)+(dx__(idx)),(y2-y1)+(dy__(idx)));
	end
end
[k,n^4]
c = c/n^4;

c
c_
end
