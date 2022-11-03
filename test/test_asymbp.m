% 2021-10-08 16:18:02.555013234 +0200

f0=10; n=1e3; x=linspace(0,1,1e3)'; y = randn(n,1);  rho = filter_f0_to_rho(f0*[1,2],1,n); y=bandpass1d_implicit(y,rho(1),1,true).^2; y(:,2) = bandpass1d_implicit(y,[rho(2),0],1,true);  y(:,3) = filter(1,[1,-rho(2)],y(:,1)); y = y./max(y); plot(y)

if (0)

n=1e4;
x=linspace(0,1,n)';
mode = 'impulse';
%mode = 'random';
mode = 'sin';
switch (mode)
case {'impulse'}
 	e = zeros(n,1);
	e(end/2) = 1;
case {'random'}
	e = randn(n,1);
case {'sin'}
	e = sin(10*2*pi*x);
end

% (1 + sin(1))^2 = 1 + 2 sin + sin(1) = 1 + 2 sin + 1 + cos(2 sin)

o = 1;
f0=10*[1,2,1/2];
rho=filter_f0_to_rho(f0,1,n);

pq = [1,0];
pq = pq/rms(pq).^2;

y=e;
o_ = 1;
for idx=1:o_
%y = 2.^o*bandpass1d_implicit(y,[rho(1),rho(1)],o,true) - 0.5.^o*bandpass1d_implicit(y,1i*[rho(2),rho(2)],o,true);
%y = o*0*filter(1-rho(1),[1,-rho(1)],y) - 0.5*
y = filter(1-rho(1),[1,-(pq(1) + pq(2)*1i)*rho(1)],y);
%y = filter(1-rho(1),[1,-rho(1)],y);
%y = y/rms(y)+1;
%y = lowpass1d_implicit(y,[0,rho(3)],o,true);
end

close all
plot(x,[e,real(y),imag(y)])

if (0)
o=4;
clf;
 e2(:,1) = ones(n,1);
 e2(n/2,1)=2;
 f0=[10,5]/2;
 rho=filter_f0_to_rho(f0,1,n);
 y(:,1) = ((sin(2*pi*f0(1)*x))+1).^o;
 y(:,2) = bandpass1d_implicit(e,[rho(1),rho(1)],o,true)+1;
y=y-min(y);
 y=y./max(y);
yy = lowpass1d_implicit(y,[0,rho(2)],1,true);
 yy=yy-min(yy);
 yy=yy./max(yy);
plot(yy);
  hold on
end

end
