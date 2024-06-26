% Tue 13 Mar 12:31:29 CET 2018

%L = [2,3];
L = [1,1];
%n = [100,200];
n = 10*[100,100];

x = L(1)/2*linspace(-1,1,n(1))';
x = x*ones(1,n(2));

y = L(2)/2*linspace(-1,1,n(2));
y = ones(n(1),1)*y;

x = flat(x);
y = flat(y);

f = {}
f{1} = @(x,y) 1 + x + 1/2*x.^2 + 1/2*x.^3
f{2} = @(x,y) 1 + y + 1/2*y.^2 + 1/2*y.^3
f{3} = @(x,y) x.*y.^2 +x.^3.*y.^2;

d2f_dxy = {}
d2f_dxy{1} = @(x,y) zeros(size(x));
d2f_dxy{2} = @(x,y) zeros(size(x));
d2f_dxy{3} = @(x,y) 2*y + 6*x.^2.*y;

Dsn = derivative_2d_mixed(n,L);

figure(1);
clf

s = 0;
for idx=1:length(f)
	val  = f{idx}(x,y);
	dval = d2f_dxy{idx}(x,y);

	dval_ = Dsn*val;

	val = reshape(val,n);
	dval = reshape(dval,n);
	dval_ = reshape(dval_,n);

	dval  = dval(s+1:end-s,s+1:end-s);
	dval_ = dval_(s+1:end-s,s+1:end-s);

	subplot(3,3,1+3*(idx-1))
	imagesc(val)
	subplot(3,3,2+3*(idx-1))
	imagesc(dval)
	subplot(3,3,3+3*(idx-1))
	imagesc(dval_)
	
	
	

	dval = flat(dval);
	dval_ = flat(dval_);

	d=dval-dval_;
	'honk'
	max(d)
	min(d)
	rms(d)
	rms(dval-dval_)

	
end




