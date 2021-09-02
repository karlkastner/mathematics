% Sun 11 Jul 16:03:45 CEST 2021
r=[-1,-1];
n=200*[1,1];
nn = prod(n);
x=randn(n);
x3=randn(2*n);
order=1;
a = [0:30:90]';

mode = 'low';
rho = 0.75;

figure(1);
clf();
for idx=1:length(a)
y = {};

% anisotropic filter, arbitrary angle
y{1} = lowpass2d_anisotropic(x,rho*[0,1],deg2rad(a(idx)),mode,order);

y{2} = lowpass2d_fft(x,rho*[0,1],deg2rad(a(idx)));

% anisotropic filter, axis-parallel
y{3} = lowpass2d_anisotropic(x2,rho*[0,1],0,mode,order);

% rotate y2
y{3} = imrotate(y{3},a(idx),'nearest','crop');

for jdx=1:length(y)
	figure(1);
	subplot(3,length(a),idx+length(a)*(jdx-1))
	f = fftshift(abs(fft2(y{jdx}))).^2;
	f_ = f;
	%f_ = lowpass2d_implicit(f,0.24);
	imagesc(f_);
	title(a(idx))
	axis equal
	axis tight

if (idx==1 && jdx==1)
	se = std(f)';
	mu = mean(f)';
	figure(2)
	plot(mu+se*[-1,0,1])
end
end % for jdx

end % for idx

