n = 1e2;
m = 10;
x = linspace(-3,-1)';
y = normcdf(x);

A  = bernstein_matrix(y,m);
Ad = bernsteind_matrix(y,m);
c  = A \ x;

if (1)
npos = 100;
z    = zeros(npos,1);
ypos = innerspace(0,1,npos)';
Apos = bernsteind_matrix(ypos,m);
c_ = lsqlin(A,x,-Apos,z);
%c = c_;
c_ = c_;
end

ni = 1e2;
yi = (0:ni-1)'/(ni-1);
Ai = bernstein_matrix(yi,m);
%xi = Ai*c;

pv = 5;
Av = vander_1d(x,pv);
cv = Av\y;


clf
subplot(2,2,1)
plot(x,y);
hold on
plot([A*c,A*c_],y,'o')
plot([Ai*c,Ai*c_],yi)
x_ = linspace(-3,3);
Av_ = vander_1d(x_,pv);
plot(x_,normcdf(x_));
plot(x_,Av_*cv);
xlim([-3,3])

Adi = bernsteind_matrix(yi,m);
subplot(2,2,2)
plot([Ad*c,Ad*c_],y);
hold on
plot([Adi*c,Adi*c_],yi);
%plot(xi,cdiff(yi)./cdiff(xi))
%plot(cdiff(xi)./cdiff(yi),yi)
%hold on
%plot(Apos*c_,ypos);

subplot(2,2,3)
xi = [Ai*c,Ai*c_];
plot(cdiff(cdiff(xi))./cdiff(yi).^2,yi)

