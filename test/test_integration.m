y=[];
f = @(y) -y;
 y0 = 2;
 k=2.^(1:10)';
 for idx=1:length(k);
 dt=1/k(idx);
 y(idx,1:6)=y0;
 for jdx=1:k(idx);
 y(idx) = (1+dt*f(y(idx,1))).*y(idx,1);
 y_ = (1+0.5*dt*f(y(idx,2)))*y(idx,2);
 y(idx,2) = (1+dt*f(y_))*y(idx,2);
 y_ = 1./(1-0.5*dt*f(y(idx,3)))*y(idx,3);
 y(idx,3) = 1./(1-dt*f(y_))*y(idx,3);
 y(idx,4) = exp(dt*f(y(idx,4)))*y(idx,4);
 y_ = exp(0.5*dt*f(y(idx,5)))*y(idx,5);
 y(idx,5) = exp(dt*f(y_))*y(idx,5);
 y_ = exp(dt*f(y(idx,6)))*y(idx,6);
 y(idx,6) = exp(dt*f(0.5*(y_+y(idx,6))))*y(idx,6);
 end;
 end;
 loglog(1./k,abs(y-y(end,:)),'.');
% log10(abs(y-0.5))
 log10(abs(y(1:end-1,:)-y(end,:)))

