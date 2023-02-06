% 2022-12-02 00:25:33.449376625 +0100
function [Shat, whp, fhp, shp, nf] = suppress_low_frequency_lobe(Shat,msk,L)
	siz = size(Shat);
	[fx_i,fy_i,frr_i] = fourier_axis_2d([1,1],siz);
	[fx,fy,frr] = fourier_axis_2d(L,siz);

	% radial density, including low-frequency components
	[Sr, fr] = periodogram_radial(Shat,L);
	Sr0 = Sr.normalized;
	% radial distribution (cdf)
	iS = cumsum(Sr.normalized)*(fr(2)-fr(1));
	iS = iS/iS(end);

	% smoothing window size ensuring convergence	
	nf = round(sqrt(find(iS>0.5,1,'first')))

	% provisional estimate of the density including the spurious low-frequency
	% components by smoothing the periodogram
	Shat_smooth = Shat;
	l_ = ceil(sqrt(numel(msk)./sum(msk(:))))+1
	% proxy for mean (4 point averge, reduced to 2 points due to symmetry)
	%Shat_smooth(1,1) = 0.5*(Shat(1,2)+Shat(2,1));
	Shat_smooth(frr_i <= l_) = max(Shat_smooth(frr_i>l_ & frr_i<=l_+1));

	% TODO gaussian smoothing window
	%Shat_smooth = circfilt2(Shat_smooth,nf);
	Shat_smooth = gaussfilt2(Shat_smooth,nf);
	%Shat_smooth = ifftshift(trifilt2(fftshift(Shat_smooth),nf));

	% re-zero mean
	Shat_smooth(1,1) = 0;

	% radial periodogram of smoothed periodogram
	Sr_smooth = periodogram_radial(Shat_smooth,L);

	% minmax of spectral energy,
	% the frequency of the minimum is used for separating the densities
	[mindx,maxdx] = minmax(Sr_smooth.normalized(2:end));
	mindx = mindx+1;	
	fhp = fr(mindx);

	% remove spurious low-frequency components by multiplication
	% with Fourier window
	if(0)
		dfrr = hypot(fx(2)-fx(1),fy(2)-fy(1));
		shp = (fhp)/sqrt(2*log(2));
		wf = 1-normpdf(frr,0,shp)/normpdf(0,0,shp);
	else
		%dfrr = hypot(fx(2)-fx(1),fy(2)-fy(1));
		shp  = (fhp)/sqrt(2*log(2));
		whp  = 1-normpdf(frr,0,shp)/normpdf(0,0,shp);	
	end

	%Shat(frr_i<=hp(idx)) = 0;
	%s = (hp(idx))/sqrt(2*log(2))
	%w = 1-normpdf(frr_i,0,s)/norm;
	Shat = whp.*Shat;

end

