% Mon  4 Oct 16:43:54 CEST 2021
function test_determine_phas_shift()
	n   = 1e4;
	m   = 10;
	lambda = 1;
	f   = 1/lambda;
	L   = 10;
	x   = innerspace(0,L,n)';
	phi = pi*x/L;
	y   = sin(2*pi*f*x + phi) + 0.25*sin(4*pi*f*x + 2*phi);
%	y   = cos(2*pi*f*x + phi) + 0.25*cos(4*pi*f*x + 2*phi);
	s = 0.1;
	[x,y,phi] = cos_random_walk(1,s,L,x(2)-x(1));
	y = y(1:end-1);
	y   = y/std(y);
	y = y-min(y);

	p = lambda/L;
	[y_mu,y_std,yy,ps] = correct_phase_shift(y,p);
	%y  = zeros(10*n,1)'
%	y = [];
%	for idx=1:m
%		y = [y; circshift(yi,10*(idx-1))];
%	end

	figure(1);
	clf
	plot(y)
	figure(2)
	clf


	ni = n/length(ps);
	yy_ = zeros(ni,length(ps));
	for idx=1:length(ps)
		id = (idx-1)*ni+(1:ni);
	subplot(2,2,1)
	plot(y(id));
	hold on
	subplot(2,2,2)
	plot(circshift(y(id),-round(ps(idx))))
	hold on		


		yy_(:,idx) = y(id); %circshift(y(id),-round(ps(idx)));
	end
	y_mu_ = mean(yy_,2);
	y_std_ = std(yy_,[],2);
%	yy_(:,idx) = y(id);


	subplot(2,2,3);
	plot(y_mu_);
	subplot(2,2,4);
	y_std = y_std/std(y_mu);
	y_mu = y_mu/std(y_mu);
	plot(y_mu);
	hold on
	plot(y((1:n/m)))
	hold on
	errorarea2(1:n/m,[y_mu-y_std,0*y_mu,y_mu+y_std],'blue')

	figure(3)
	subplot(2,2,1)
	fx = fourier_axis(x(1:ni));
	plot(fx,periodogram(y_mu-mean(y_mu),L))
	xlim([0,5*f])
end

