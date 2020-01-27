% Mon  9 Dec 19:05:11 +08 2019
function ys = mysmooth(y,d,wrap)
fmode = 'splines';
% generic hydrograph by smoothing / filtering
%	zs_shifted = circshift(zs,dd)
%	plot(zs_shifted)
%	zs_shifted = reshape(zs_shifted,365,[]);
%	datetick();

switch(fmode)
case {'geo'}
	m = nangeomean(zs_shifted,2);
	r = log(zs_shifted)-log(m);
	r = m.*nanstd(r,0,2);
case {'lin'}
	m = nanmean(zs_shifted,2);
	s = nanstd(zs_shifted,[],2);
case {'linfft'}
	% TODO, there is a sinosoidal seasonal component and a recession component
	%F    = fourier_matrix(365*[1,1/2,1/4],1:365);
	%-> here lsqnonline with filter for recession
	%-> one leading year as buffer
	nf   = 8;
	F    = fourier_matrix(365*1./(1:nf),1:365);
	c    = F \ h_p;
	cQ   = F \ Q_p;
	c_ = fourier_cesaro_correction(c,nf);
	cQ_ = fourier_cesaro_correction(cQ,nf);
	if (0) % fourier with recession
		fun = @(c) filter(1-c(end),[1,-c(end)],F*c(1:end-1),h_p(1))
		c(end+1) = 0.999;
		c = lsqnonlin(@(c) fun(c) - h_p,c);
	end
case {'ptn'}
	r = 0.95;	
	chf = [3e4, 1, 2.4, 2e3, r, 0*3*r^2, 0*-r^3, r, 0*3*r^2, 0*-r^3];
	chf = [3e4/1e4, 1.148135  2.370936, 2e3/1e3, 0.5258    0.1859    0.1807    0.5403    0.1800    0.1833]
	chf = [2.7932   1.8578    1.2125    2.3690    0.4736    0.1661    0.2570    0.4163    0.2780    0.2100];
	chf = lsqnonlin(@(c) hfilter(c(1),c(2),c(3),c(4),c(5:7),c(8:10),365,Q_p) - Q_p, chf);
case {'hfilter'}
	h_p__ = hfilter(chf(1),chf(2),chf(3),chf(4),chf(5:7),chf(8:10),365,Q_p);
case {'splines'}
	%	[h_p4,a,b] = fit(cvec(1:365),Q_p(:,2),'smoothingspline','SmoothingParam',0.07);
	ys = smooth_with_splines(y,d,wrap);
if (0)
	r2 = (zs(1:nd)-zs_p).^2.*(ny)/(ny-1);
	r2 = sqrt(r2);
	%r2 = r2./zs_p.^2;
	F=fourier_matrix(365*[1,1/2],1:nd);
	cr2 = F\r2;
	r2p = (F*cr2);
	dc = min(0,min(r2p))
	cr2(1) = cr2(1)-dc+sqrt(eps);
	r2p = (F*cr2);
	%r(end+1:365*ny) = NaN;
	%r = nanstd(reshape(r,365,[]),[],2);
	m = nanmean(zs_shifted,2);
end
end % switch
	% TODO residuals
end 

