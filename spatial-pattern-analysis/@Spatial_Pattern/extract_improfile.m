% alternative extraction to rotating the image

	s  = 4;
	ic0 = [n(1)/2+1, imax];
	jc0 = [n(2)/2+1, jmax];
	ic = ic0(1) + [0, s*(ic0(2)-ic0(1))]-1*0;
	jc = jc0(1) + [0, s*(jc0(2)-jc0(1))]-1*0;

	[xc_,yc_,cS] = improfile(S,round(jc),round(ic));
	[xc_,yc_,cR] = improfile(R,round(jc),round(ic));

	n = size(S);  % was R'
	r = hypot(pg.x,rvec(pg.y));
	rc  = r(sub2ind(n,round(xc_),round(yc_)));
	%fr  = fourier_axis(rc);
	frc = fr(sub2ind(n,round(xc_),round(yc_)));
	
	% central frequency
	[~, mdx] = max(cS);
%	fc(kdx) = frc(mdx);

	% perpendicular density
	ic0_ = [n(1)/2+1, imax];
	jc0_ = [n(2)/2+1, jmax];
	ic_ = 0*ic(1) + imax + [0, -s*(jc0_(2)-jc0_(1))];
	jc_ = 0*jc(1) + jmax + [0,  s*(ic0_(2)-ic0_(1))];

	% cross density
	[xc_,yc_,aS] = improfile(S,round(jc_),round(ic_));
	fxx = repmat(cvec(fx),1,length(fy));
	fyy = repmat(rvec(fy),length(fx),1);
	fx_ = fxx(sub2ind(n,round(xc_),round(yc_)));
	fy_ = fyy(sub2ind(n,round(xc_),round(yc_)));
	fa = hypot(fx_-fx_(1),fy_-fy_(1));
	aS = aS/(sum(aS)*(fa(2)-fa(1)));
%diff(fx(round(ic_)))*diff(fx(round(ic))) + diff(fy(round(jc_)))*diff(fy(round(jc)))

	xc_ = [n(1)/2+1, pg.stat.imax];
	yc_ = [n(2)/2+1, pg.stat.jmax];
	d   = [xc_(2)-xc_(1);
	       yc_(2)-yc_(1)];
	
	d = d./hypot(d(1),d(2));
	if (0)
	xc_ = [1,n(1)];
	yc_ = xc_.*d(2)/d(1);
	else
		yc_ = [1,n(2)];
		xc_ = -fliplr(yc_)*d(1)/d(2);
		xc_ = xc_-mean(xc_) + n(1)/2;
	end
	
	
	% transect and 1D estimates along transects
	[xc_,yc_,cI] = improfile(pg.b',round(xc_),round(yc_));
	fdx_ = isfinite(cI);
	xc_ = xc_(fdx_);
	yc_ = yc_(fdx_);
	cI = cI(fdx_);
	n1D = length(cI);
%plot([xc_,yc_,cI])
%pause
	cI = cI-mean(cI);
%	cI = cI/rms(cI);
	rI   = hypot(xc_-xc_(1),yc_-yc_(1));
	L1D  = rI(end);
	df1D = 1./L1D;
	f1D  = fourier_axis(rI);
	%x1D = fourier_axis(rI);
	Shat1D = abs(fft(cI)).^2;
	% normalize
	Shat1D = 2*Shat1D/(sum(Shat1D)*df1D);

	% autocorrelation
	%R1D    = autocorr_fft(cI);
	R1D    = 0.5*n1D/L1D*real(ifft(Shat1D));

	% smoothing window	
	x1D = linspace(-L1D/2,L1D/2,n1D)'; % = rI - L1D/2
	w1D = radial_window(x1D,rmax);

	%S1D = periodogram_bartlett(cI,L1d,nb,n1d)
	%S1D(:,2) = S1D(:,1).^2;
	S1D = 2*L1D/n1D*real(fft(ifftshift(w1D).*R1D));

	
disp('test')
size(cI)
sum(flat(Shat)*df(1)*df(2))
sum(flat(S)*df(1)*df(2))
max(R(:))
sum(Shat1D*df1D)
sum(S1D*df1D)
sum(sqrt(cS)*df1D)
R1D(1)
L1D*pg.stat.fcr

if (0)
figure(1)
clf
subplot(2,2,1)
imagesc(pg.b)
subplot(2,2,2)
%img = img-min(img);
%img = img./max(img);
%p = 0.5;
%imagesc((1-(1-img).^p).^p)
colormap gray
subplot(2,2,3)
imagesc(pg.b)
caxis(quantile(pg.b(:),[0.1,0.7]))
end

