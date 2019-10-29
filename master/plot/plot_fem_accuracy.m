% Wed Aug 29 01:09:44 MSK 2012
% Karl Kästner, Berlin

	lw  = 1;
	lw_ = 1;
	ms = 12;

	a = 0;
	b = 2;
	y_max = ceil(exp(b));
	X_true = linspace(a,b,100)';
	Y_true = exp(X_true);


	subplot(3,3,1); cla();
	X = a + (b-a)*[0 1];
	Y = exp(X);
	plot(X_true,Y_true,'-r','linewidth',lw,'markersize',ms); hold on
	plot(X,Y,'.-b','linewidth',lw_,'markersize',ms);
	axis([a b 0 y_max]);
	axis square
	xlabel('x'); ylabel('y');
	p=get(gca,'position')
	set(gca,'position',[p(1) 0.75*p(2) p(3) p(4)])
	title([ 'h=2 p=1' ] );

	subplot(3,3,2); cla();
	X = a + (b-a)*[0 0.5 1];
	Y = exp(X);
	plot(X_true,Y_true,'-r','linewidth',lw,'markersize',ms); hold on
	plot(X,Y,'.-b','linewidth',lw_,'markersize',ms);
	axis([a b 0 y_max]);
	axis square
	xlabel('x'); ylabel('y');
	title([ 'h=1 p=1' ] );

	subplot(3,3,3); cla();
	X = a + (b-a)*[0 0.25 0.5 0.75 1];
	Y = exp(X);
	plot(X_true,Y_true,'-r','linewidth',lw,'markersize',ms); hold on
	plot(X,Y,'.-b','linewidth',lw_,'markersize',ms);
	axis([a b 0 y_max]);
	axis square
	xlabel('x'); ylabel('y');
	title([ 'h=1/2 p=1' ] );

	subplot(3,3,5); cla();
	X = a + (b-a)*[0 0.5 1]';
	Y = exp(X);
	A = [X.^0 X X.^2];
	c = A \ Y;
	Y_ = [X_true.^0 X_true X_true.^2]*c;
	plot(X_true,Y_true,'-r','linewidth',lw,'markersize',ms); hold on
	plot(X,Y,'.b','linewidth',lw,'markersize',ms);
	plot(X_true,Y_,'-b','linewidth',lw_);
	axis([a b 0 y_max]);
	axis square
	xlabel('x'); ylabel('y');
	title([ 'h=1 p=2' ] );

	subplot(3,3,6); cla();
	X = a + (b-a)*[0 1/3 2/3 1]';
	Y = exp(X);
	A = [X.^0 X X.^2 X.^3];
	c = A \ Y;
	Y_ = [X_true.^0 X_true X_true.^2 X_true.^3]*c;
	plot(X_true,Y_true,'-r','linewidth',lw,'markersize',ms); hold on
	plot(X,Y,'.b','linewidth',lw,'markersize',ms);
	plot(X_true,Y_,'-b','linewidth',lw_);
	axis([a b 0 y_max]);
	axis square
	xlabel('x'); ylabel('y');
	title([ 'h=1 p=3' ] );

	print -depsc ../img/fem-approx.eps
	system('epstopdf ../img/fem-approx.eps')

