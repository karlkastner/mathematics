% 2013-07-25 23:28:43.000000000 +0200

%compare to fft, different step width, formulate

% 1st : crosstalk, especially to higher frequencies
% 2nd : frequency depending scaling issue	-> inverse does not reproduce exact amplitude
% 3rd : cone of inlfuence at the boundary

% TODO, third frequency and different amplitudes
% TODO, reduce to cone of influence
% TODO, plot time series
% TODO, manual transform, convolve sine with gausian

function test_wavelet()
	omega = [1 0.5 0.25];
	A     = [1 1 1];
	phi   = [0 0 0];
	
	% todo, vary phi2 in loop from 0 to 2pi

	n = 100;
	k = 20; %*max(omega)/min(omega);
	t = linspace(0,k*2*pi,n)';

	x(:,1) = harmonics(t,A(1),omega(1),phi(1));
	x(:,2) = harmonics(t,A(2),omega(2),phi(2));
	x(:,3) = harmonics(t,A(1:3),omega(1:3),phi(1:3));
	%x12= harmonics(t,A(1:2),omega(1:2),phi(1:2));
	x(:,4) = harmonics(t,A(1:2),omega(1:2),[0 pi/2]);

	scales = 1:0.5:n/4;
	f = [];
	wname = 'morl';
%wname = 'haar';
	for idx=1:4
		c(:,:,idx) = cwt(x(:,idx),scales,wname);
		f(:,idx) = fft(x(:,idx));
	end
	c1 = c(:,:,1);
	c2 = c(:,:,2);
	c12 = c(:,:,3);
	c12_ = c(:,:,4);
	
	m1 = max(abs(c1),[],2);
	m2 = max(abs(c2),[],2);
	m12 = max(abs(c12),[],2);
	m12_ = max(abs(c12_),[],2);

	subplot(3,2,1);
	cla();
	plot(scales,[m1 m2 m12 m12_]);
	hold on;
	plot(10*abs(f(1:end/2,1:2))/max(abs(f(:))),'--');

	m1   = sqrt(sum(c1.^2,2)/n);
	m2   = sqrt(sum(c2.^2,2)/n);
	m12  = sqrt(sum(c12.^2,2)/n);
	m12_ = sqrt(sum(c12_.^2,2)/n);

	subplot(3,2,3)
	plot(scales,[m1 m2 m12 m12_]);


	m1   = sqrt(sum(abs(c1),2)/n);
	m2   = sqrt(sum(abs(c2),2)/n);
	m12  = sqrt(sum(abs(c12),2)/n);
	m12_ = sqrt(sum(abs(c12_),2)/n);

	subplot(3,2,5)
	plot(scales,[m1 m2 m12 m12_])

	tick = get(gca,'xtick');
%	t = 2+pi*2.^(0:10);
%	set(gca,'xtick',t);
	l = 1;
	l = 2*pi;
	set(gca,'xticklabel',l./(tick))

	M = [5 10 20];
	for idx=1:3
		y = zeros(size(x,1),1);
		for jdx=1:3
			%m = round(2*pi/omega(jdx))
			%m = M(jdx);
%			m = round(n/(k)*1/(omega(jdx)))
			m = round(n/(2*pi)*omega(jdx))
			%y = [y c(m,:,idx)'];
			y = y + c(m,:,idx)';
		end
		subplot(3,2,2*idx);
		cla();
		plot(y);
		hold on
		plot(x(:,idx),'--');
	end

end

function x = harmonics(t,A,omega,phi)
	x = zeros(size(t));
	for idx=1:length(A)
		%x = x + A(idx)*sin(omega(idx)*t + phi(idx));
		x = x + A(idx)*sign(sin(omega(idx)*t + phi(idx)));
	end
end


