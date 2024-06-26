% Thu  1 Sep 12:13:52 CEST 2016

function test_line_search_quadratic2()

%x0 = 2;
%w = ([1 10])';
w = (1:10)';
x0 = 10*ones(length(w),1);

opt.lsfun = @line_search_quadratic2;
opt.lsmaxiter = 1000;
opt.verbose = 0;
opt.lsh = 0.1;
opt_ = opt;
opt_.lsfun = @line_search_wolfe;
%opt__ = opt;
%opt__.lsh = 1e-3;
%opt__.lsfun = @line_search_quadratic;
%order = 3;
%opt___ = opt;
%opt___.lsfun = @(fun,x,f0,g,dir,h,lb,ub,maxiter,v) line_search_polynomial(fun,x,f0,g,dir,h,lb,ub,maxiter,order,v)


%P = (1.1:0.1:10);
P = linspace(0,1,100);
%P = (1.2:0.01:4);
N = [];
M = [];
X = [];
if (1)
for idx=1:length(P)
idx
p = P(idx)
cnt = 0;
[X(idx,:) v1 v2 M(idx,1)] = nlcg(@fun,x0,opt);
N(idx,1) = cnt;
cnt = 0;
[Y(idx,:) v1 v2 M(idx,2)] = nlcg(@fun,x0,opt_);
N(idx,2) = cnt;
%cnt = 0;
%[Z(idx,:) v1 v2 M(idx,3)] = nlcg(@fun,x0,opt__);
%N(idx,3) = cnt;
%cnt = 0;
%[Z(idx,:) v1 v2 M(idx,4)] = nlcg(@fun,x0,opt___);
%N(idx,4) = cnt;
end

subplot(2,2,1)
plot(P,N)
legend('Quadratic','Wolfe','local quadratic')
hline(median(N(:,1)))
hline(median(N(:,2)),'r')
subplot(2,2,2)
plot(P,M)
subplot(2,2,3)
plot(P,X)
end

cnt = 0;
p = 0.55;
nlcg(@fun,x0,opt); 

x = linspace(-5,5,100)+1;
y = linspace(-5,5,100)+1;
subplot(2,2,4)
z=[];
for idx=1:length(x)
for jdx=1:length(y)
%	z(idx,jdx) = fun([x(idx); y(jdx)]);
end
end
%imagesc(x,y,z)
%plot(x,fun(x))


function [f g] = fun(p,x)
	cnt = cnt+1;
%	f = w'*(x-1).^2;
	f = w'*(p*abs(x-1).^4) + w'*(1-p)*abs(x-1).^2;
%	f = w'*(p*abs(x-1).^4) + (max(w)+1-w)'*(1-p)*abs(x-1).^2;
%	f = w'*(abs(x-1).^p);
%	f = f + w'*((x-1).^2);
	if (nargout() > 1)
		g = grad(@fun,x,[],'one-sided');
	end
end

end
