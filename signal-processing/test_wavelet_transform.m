% Fri  8 Jul 16:41:06 CEST 2016
% Karl Kastner, Berlin

dt = 1/(3*24);
T  = 100;
t  = cvec(0:dt:T);

f0 = 1;
n = 7;

x  = 1*sin(2*pi*f0*t);
y = wavelet_transform_single(x,dt,f0,n,'flattopwin');

figure(1);
subplot(2,2,1);
plot(t,[x,abs(y)]);
subplot(2,2,2);
plot(t,[angle(y)]);

% frequency response
f = f0*logspace(-1,1,1000)';
A = zeros(size(f));
winstr = {'boxcar','triangular','hanning','flattopwin'};
for idx=1:length(f)
	x  = 1*sin(2*pi*f(idx)*t);
	for jdx=1:length(winstr)
		y = wavelet_transform_single(x,dt,f0,n,winstr{jdx});
		A(idx,jdx) = abs(y(round(end/2)));
	end
%	figure(3)
%	subplot(3,4,idx);
%	plot(abs(y))
end
figure(2)
loglog(f,A)
for idx=-2:2
	vline(f0*2.^idx);
end
legend(winstr);

