% 2022-12-02 00:25:33.449376625 +0100
function obj = analyze_grid_2(img,msk,ns)
	area_msk_pxl2 = sum(msk(:));
	% apply mask
	mu = mean(img(msk));
	img = img - mu;
	img(~msk) = 0;


	% original image before making square
	b     = img;
	bmsk  = msk;

	siz = size(img);
	% make square
 	if (siz(1) < siz(2))
		siz(1) = siz(2);
		img(siz(1),1) = 0;
		msk(siz(1),1) = false;
	end
	if (siz(1) > siz(2))
		siz(2) = siz(1);
		img(1,siz(2)) = 0;
		msk(1,siz(2)) = 0;
	end
	% dummy dimension for images, as it is for most images not known
	% since we are only interested in regularity, the actual length-scale is irrelevant
	% note, since the image has been made square, this is [1,1]
	%L = siz./geomean(siz);
	L = [1,1];
	[fx,fy,frr] = fourier_axis_2d(L,siz);

	% 2D periodogram (not normalize)
	Shat = abs(fft2(img)).^2;
	
	% removal of spurious low frequency components
	% we assume that spurious low frequency components are approximately isotropic
	[Shat, wlp, lp] = suppress_low_frequency_lobe(Shat,msk,L);

	% radial density (radial periodogram)
	[Sr, fr] = periodogram_radial(Shat,L);
	Sr = Sr.normalized;

	% cumulative distribution
	iSr = cumsum(Sr);
	iSr = iSr/iSr(end);
	% make values unique (quick hack)
	iSr = cvec(iSr) + (0:length(iSr)-1)'*1e-12; 
	f_05 = interp1(iSr,fr,0.05,'linear');
	f_50 = interp1(iSr,fr,0.50,'linear');
	f_95 = interp1(iSr,fr,0.95,'linear');

	dfr = frr(1,2);
	nf   = round(sqrt(f_05/dfr));
%	nf   = round(sqrt(find(iSr>0.5,1,'first')))
%	nf   = round(sqrt(find(iSr>0.5,1,'first')));
	Sbar = gaussfilt2(Shat,nf);


	dfrr = sqrt(frr(2,1)*frr(1,2));
	% upper frequency range limit for the periodicity test
	%frmax = f_95+dfrr;
	% lower frequency range lomot for the periodicity test
	%frmin = max(f_05-dfrr,lp*dfrr)
	%frmin = max(0.2*frmu,lp(idx)*dfrr);
	df    = sqrt(L(1)*L(2));
	%nf  = round(0.5*frmu/df);
	%nft = round(0.25*(f_95-f_05)/dfrr);
	% smoothing window radius for frequency test
	nft  = round(0.25*f_50/dfrr);
	nft  = max(nft,2);

	% repeat for non-square image b
	if (0)
	% TODO as the mask is used, it is not necessary to repeat it with a square image
	sizb = size(b);
	Lb = sizb / sqrt(sizb(1)*sizb(2));
	[fbx,fby,fbr] = fourier_axis_2d(Lb,size(b));
	Shat_b = abs(fft2(b)).^2;
	% TODO anisotropic smoothing due to non-square image
	Sbar_b = gaussfilt2(Shat_b,nf);
	%Sbar_b = ifftshift(trifilt2(fftshift(Shat_b),nf));
	Sbar_b(fbr < log(4)*slp) = 0;
	else
		b    = img;
		sizb = size(b);
		Lb = sizb / sqrt(sizb(1)*sizb(2));
		bmsk = msk;
		Sbar_b = Sbar;
		fbx = fx;
	end
		%Sbar_b(fbr < log(4)*slp) = 0;
		[Ssort,sdx] = sort(Sbar_b(:));
		iSort = cumsum(Ssort);
		% restrict test to frequency range that contains upper 90% of spectral energy
		iSsort = cumsum(Ssort);
		fdx = (iSsort >= 0.2*iSsort(end));
 		fmsk = false(size(b));
		fmsk(sdx(fdx)) = true;
		fmsk = fmsk & (frr > dfr*lp);
		% by symmetry, excluded the symmetric half plane
		fmsk = fmsk & cvec(fbx) >= 0;
if (0)
figure(1e3)
clf()
lp
subplot(2,2,1)
imagesc(Shat)
subplot(2,2,2)
imagesc(Sbar)
subplot(2,2,3)
imagesc(fmsk)
contour(fmsk,[-1,0.5])
colorbar
pause
end
	if (numel(bmsk)*ns<2e8)
		try
		%[p_periodic,stati,ratio,qq] = periodogram_test_periodicity_2d(b, Lb, nft, bmsk, fmsk, ns);
		%[p_periodic,stati,ratio,qq] = periodogram_test_periodicity_2d(b, Lb, nft, bmsk, frmin, frmax, ns);
		[p_periodic,stati,ratio,qq] = periodogram_test_periodicity_2d(b, Lb, nft, bmsk, fmsk, ns);
		%[p_periodic,stati,ratio,qq] = periodogram_test_periodicity_2d(b, Lb, nft, bmsk, fmsk, ns);
		fr_periodic = stati.fr_max;
		catch e
			e
			p_periodic = NaN;
			fr_periodic = NaN;
		end
	else
		p_periodic = NaN;
		fr_periodic = NaN;
	end
	%nf_a(idx,1) = nf;
	%frmu = f_50;	

	[S_,f_,piso]=separate_isotropic_from_anisotropic_density(Shat,L,nf);
	isisotropic = piso>=0;

	Sx = mean(S_.hat,2);
	Sx(f_.x<0) = 0;
	Sx = Sx/(sum(Sx(f_.x>=0)));
	Sy = mean(S_.hat,1);
	Sy(f_.y<0) = 0;
	Sy = Sy/(sum(Sy(f_.y>=0)));

%	if (0)
%		nf = round( 1/3*(find(iS>0.75,1,'first')-find(iS>0.25,1,'first')));
%	else
		% convergent criterium for L->inf and dx->0
		% reduced smoothing bandwidth
		% was 1/3
%		nf_ = round(1/2*(find(iS>0.75,1,'first')-find(iS>0.25,1,'first')));
%		if (nf_<nf)
%			nf = nf_;
%		end
%	end	


	if (0) %nf > 1)
		Sr(2:end) = trifilt1(Sr(2:end),nf);
		Sx(2:end) = trifilt1(Sx(2:end),nf);
	end

if (0)
	% automati.pdfy determine highpass filter length, skip zero mean component
	if (isisotropic)
		kri = cvec(0:(length(Sr)-1));
		[mindx,maxdx] = minmax(kri(2:end).*Sr(2:end));
	else
		[mindx,maxdx] = minmax(Sx(2:end));
	end
	lp_ = mindx+1;

	% remove spurious low-frequency compoents by highpass filtering
	Shat(frr_i<lp_)       = 0;
	Sx(1:lp_)  = 0
	Sr_ = Sr;
	Sr(1:lp_)      = 0;
	% renomalize
	dfr = fr(2)-fr(1);
	Sr = Sr./(sum(Sr)*dfr);
	dfx = fx(2)-fx(1);
	Sx = Sx./(sum(Sx)*dfx);
end % if 0
	
	if (isisotropic)
		[Sc,ldx] = max(Sr);
		lc       = 1./fr(ldx);
	else
		[Sc,ldx] = max(Sx);
		lc       = 1./abs(fx(ldx));
	end
	
	regularity = Sc./lc;

	obj = struct();
	obj.area_msk_pxl2 = area_msk_pxl2;
	obj.fr_periodic  = fr_periodic;
	obj.f_50         = f_50;
	obj.isisotropic  = isisotropic;
	obj.p_isotropic = piso;
	obj.lc       = lc;
	obj.lp       = lp;
	obj.nf       = nf;
	obj.nft     = nft;
	obj.p_periodic   = p_periodic;
	obj.regularity   = regularity;
	obj.Sc           = Sc;
	obj.siz      = siz;

	obj.Sr = Sr

	obj.fr_05 = f_05;
	obj.fr_50 = f_50;
	obj.fr_95 = f_95;


	obj.Sr = Sr;
	obj.Sx = Sx;
	obj.fr = fr;
	obj.fx = fx;
	obj.fy = fy;
	obj.Sx = Sx;
	obj.fmsk = fmsk;
	obj.Srot = S_;
end % analyze_grid2

