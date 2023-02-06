% 2023-01-15 22:45:58.576936832 +0100

	[Rr,r,Rt,t,Rc,xc,yc,out] = autocorrelation_max(R,L)

	n = size(R);
	clf

	% radial
	subplot(2,2,1);
	imagesc(R);
	hold on;
	plot(out.rij(2,:),out.rij(1,:),'r')
	plot(out.jc+floor(n(2)/2)+1,ic+floor(n(1)/2)+1,'o');

	subplot(2,2,2)
	plot(Rr);

	% angular
	subplot(2,2,3)
	imagesc(R)
	% max
	plot3(out.jc,out.ic,1,'r*')
	% start point
	plot(out.tij(2,1),out.tij(1,1),'r*');
	% ring
	plot(out.tij(2,:),out.tij(1,:),'r*');

	subplot(2,2,4)
	plot(Rt)


