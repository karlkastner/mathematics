% 2013-07-26 01:51:46.000000000 +0200
% Karl Kastner, Berlin
%
% wavelet transform (parseval, scaling)
% inverse wavelet transform
% relation to fourier transform
%
function test_wavelet()
	omega = [1 0.5 0.25];
	A     = [0 1 0];
	phi   = [0 0 0];
	
	% todo, vary phi2 in loop from 0 to 2pi

	n = 1e3-1;
	k = 10; %*max(omega)/min(omega);
	t = linspace(0,k*2*pi,n)';
	x = zeros(size(t));

	for idx=1:length(A)
		x = x + A(idx)*sin(omega(idx)*t + phi(idx));
	end

%	figure(1); clf
%	subplot(2,2,1)
%	plot(t,x)
%	c = cwt(A,1:n,'morl');
%	subplot(2,2,2)
%	imagesc(abs(c)')
%	colormap gray;

% problems at the boundary as long as one periode

%	scales = 2.^(-5:5);
%	scales = (1:n);
%	scales = 1:2:n;
%	scales = 1:2*n; %0.1:0.1:1000; %2*pi*omega;
%	scales = logspace(0,log(n)/log(10),100);
%	scales = n*omega;
%	scales  = n*[0.125 0.25 0.5 1 2]
%	wname  = 'db1';	% what is the number for ?
%	wname  = 'gaus';
	c = [];
%	[ca cd] = dwt(x,'haar');
%	plot([x [ca; cd]]);
%	for idx=1:9
%		subplot(3,3,idx)
%		[x cd] = dwt(x,'haar');
%		y(idx) = norm(cd);
%		plot(x)
%		wname = ['haar' num2str(idx)]
%		c = [c dwt(x,wname)];
%	end
%	plot(y)

	scales = 1:n;
	wname  = 'haar';
	wname  = 'morl';
	[c ] = cwt(x,scales,wname); %,'scal');
%	plot(c);
	figure(1); clf
	subplot(2,2,1)
	plot(t,x)
	xlim([t(1) t(end)])
	ylim([min(x) max(x)])
	subplot(2,3,4)
	imagesc(t,scales,abs(c))
	subplot(2,3,5)
	imagesc(t,scales,atan2(real(c),imag(c)))
	%imagesc(abs(log(c)))
	colormap gray
	subplot(2,3,6)
	plot([max(abs(c),[],2) max(real(c),[],2)])
%	plot(sum(abs(c),2)/size(c,2))
pause
	% reconstruct
	y = zeros(size(x));
	A   = abs(c);
	phi = atan2(imag(c),real(c)); % + rad2deg(180);
%	phi = atan2(real(c),imag(c));
	[v mdx] = max(A(:));
%	A = A(2:end-1:2:end-1);
	[idx jdx] = ind2sub(size(A),mdx)
	length(scales)/k*scales(idx)

	ss = scales(idx)
	sc = length(scales)/(k*scales(idx))
%	ss = scales(idx)/(t(end) - t(1))
	for idx=1:length(scales)
		for jdx=1:length(t)
			y(jdx) = y(jdx) + A(idx,jdx)*sin( length(scales)/(scales(idx)*k)*t(jdx) + phi(idx,jdx));
%			y(jdx) = y(jdx) + A(idx,jdx)*sin( scales(idx)*t(jdx) + phi(idx,jdx));
%			y(jdx) = y(jdx) + A(idx,jdx)*sin( (n/k)*1/scales(idx)*t(jdx) + phi(idx,jdx));
		end
	end
	norm(y)
	norm(A)
	y = y/norm(y)*sqrt(length(y));
	subplot(2,2,1);
	hold on;
	size(y)
	plot(t,y,'g');
	hold off
	axis auto


end % function test_wavelet()

