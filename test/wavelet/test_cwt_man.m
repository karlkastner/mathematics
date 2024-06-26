% 2013-07-26 01:51:11.000000000 +0200
% Karl Kastner, Berlin
%
	
	n = 1e3-1;
	k = 2*pi*16;
	f = 1;
	omega = 2*pi*f;
	phi = 0;
	t = linspace(0,k,n)';
	A = 1;

	figure(1);
	clf();
	K = -4:1;
	for kdx=1:length(K)+1
		col = 0.5*(1+[ sin(2*pi*(kdx-1)/(length(K)+1)) sin(2*pi*(kdx-1)/(length(K)+1)+pi/3) sin(2*pi*(kdx-1)/(length(K)+1) +2*pi/3)]);


%	x = sin(omega*t + phi) + sin(2*omega*t + phi) + sin(4*omega*t + phi) + sin(8*omega*t) +s
	if (kdx < length(K))
		x = A*sin(2^K(kdx)*omega*t + phi);
	else
		x = zeros(size(t));
		for idx=1:length(K)
			x = x + A*sin(2^K(idx)*omega*t + phi);
		end
	end
%	x = randn(n,1)-0.5;
		
	scales = 1:n;
	scales = logspace(0,log10(n),100);
	c = cwt(x, scales, 'morl');
	scales_ = scales;
	c = c .* sqrt((1./scales_)'*ones(size(x))');
	sigma = 5;
	c_man = cwt_man(x,scales,sigma);
	c_man = c_man .* sqrt((1./scales_)'*ones(size(x))');
%	s = 0.1:0.1:20;
%	for sdx=1:length(s);
%		c_man = cwt_man(x,s(sdx)*scales);
%		y(sdx) = norm(c-c_man);
%	end
%	c = c/max(c(:));
%	c_man = c_man/max(c_man(:));

%	c = c/sqrt(n);
%	c_man = 2*c_man / n;

%	figure(2)
%	plot(s,y);

%	figure(1);
%	clf();
	subplot(3,2,1);
	imagesc(c);
	colorbar()
	caxis([0 3]);
	colormap gray;
	subplot(3,2,3);
	imagesc(c_man);
	caxis([0 3]);
	colorbar()
	subplot(3,2,5);
%	imagesc(((c./c_man)));
	imagesc(c - c_man);
	colorbar()
	caxis([-1.5 1.5]);

	subplot(3,2,2);
%	plot(sqrt(sum(c.^2,2)/n),'color',col); hold on
	plot(max(c,[],2),'color',col);
	set(gca,'xscale','log')
	%plot(sum(abs(c)/n,2),'color',col); hold on
	%plot(max(c,[],2))
	ylim([0 1.5])
	subplot(3,2,4);
%	plot(sqrt(sum(c_man.^2,2)/n),'color',col); hold on
	plot(max(c_man,[],2),'color',col);
	ylim([0 1.5])
	set(gca,'xscale','log')
	%plot(max(c_man,[],2))
%	caxis(1e-1*[-2 2])

end % for kdx

