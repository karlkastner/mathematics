% Wed 18 May 13:50:47 CEST 2022
% 2022-12-02 00:25:33.449376625 +0100
function obj = analyze_grid(obj)
	% output
	S    = obj.S;
	R    = obj.R;
	stat = obj.stat;

	n = obj.n;
	L = obj.L;
		
	b   = obj.b;
	bmsk = obj.msk.b;
	if (isempty(bmsk))
		bmsk = true(size(b));
	end

	% convert image to double
	b    = double(b);
	% convert mask to logical
	bmsk = (bmsk > 0);

	% fraction of ground coverged
	% there is a matlab bug that double images need to be scaled
	bscaled = (b - min(b(bmsk)))/(max(b(bmsk))-min(b(bmsk)));

	thresh_b      = graythresh(bscaled(bmsk));
	stat.coverage = sum(bscaled(bmsk)<thresh_b)/sum(bmsk(:));

	% resample, to make dx identical to dy
	% (depending on the grid projection, this might not be the case)
	dx = L./n;
	n_ = round(L./min(dx));
	if (n_(1)>n(1))
		% TODO implement
		% note : this is not necessary as the current input files have dx=dy	
		error('not yet implemented');
	end
	if (n_(2)>n(2))
		% TODO implement
		% note : this is not necessary as the current input files have dx=dy	
		error('not yet implemented');
	end
	n = n_;

	% resize, to make domain square
	nmax = max(n);
	L    = L.*nmax./n;
 	if (n(1) < nmax)
		n(1) = nmax;
		b(nmax,1) = 0;
		bmsk(nmax,1) = false;
	end
	if (n(1) > n(2))
		nmax = n(1);
		b(1,nmax) = 0;
		bmsk(1,nmax) = 0;
	end
	df   = 1./L;
	n    = [nmax,nmax];
	dx   = L./nmax;
	stat.L_square = L;
	stat.n_square = n;
	obj.b_square = b;

	% grid in real space
	% note : x an y do not match original figure when the figure was resampled or or resized
	obj.x   = linspace(-L(1)/2,L(1)/2,n(1))';
	obj.y   = linspace(-L(2)/2,L(2)/2,n(2))';
	%rr      = hypot(obj.x',obj.y);

	% grid in frequency space
	[obj.f.x,obj.f.y,frr,ftt] = fourier_axis_2d(L,n);

	% statistics of the masked area
	area_msk = sum(bmsk(:))*dx(1)*dx(2);
	% extend approximation by an ellipsis
	nmsk = sum(bmsk(:));
	centroid.x = sum(bmsk*cvec(obj.x))./nmsk;
	centroid.y = sum(rvec(obj.y)*bmsk)./nmsk;
	sx2 = sum(bmsk*(cvec(obj.x)-centroid.x).^2)./nmsk;
	sy2 = sum(bmsk*(cvec(obj.y)-centroid.y).^2)./nmsk;
	sxy = sum(bmsk*((cvec(obj.x)-centroid.x).*(cvec(obj.y)-centroid.y)))./nmsk;
	centroid.C = [sx2,sxy;
        	      sxy,sy2];
%	centroid.Crot = rotmat(-angle_deg)*centroid.C;
%	centroid.L = sqrt([centroid.Crot(1,1),centroid.Crot(2,2)]); 
% anti rotation for area
%r = sxy./sqrt(sx2*sy2)
%Rot = [1,sxy;
%       -sxy,1];
%centroid.sr = sqrt(sx2 + 2*sxy + sy2);


	% the weighted mean allows for feathering the transition
	mu  = wmean(double(bmsk(:)),b(:));

	% subtract mean
	b  = (b - mu);

	% normalize
	b = b/wrms(bmsk(:),b(:));

	% 2D periodogram, not yet normalized
	S.hat  = abs(fft2(bmsk.*b)).^2;

	% removal of spurious low frequency components
	% we assume that spurious low frequency components are predominantly isotropic
	[S.clip, ~, fhp, shp] = suppress_low_frequency_lobe(S.hat,bmsk,L);
	%[S.clip, w] = suppress_low_frequency_components_1(obj,S.hat)


	% radial density (radial periodogram)
	[Sr, obj.f.r] = periodogram_radial(S.hat,L);
	S.r = Sr.normalized;

	obj.w.x = 1-normpdf(obj.f.x,0,shp)/normpdf(0,0,shp);	
	obj.w.y = 1-normpdf(obj.f.y,0,shp)/normpdf(0,0,shp);	
	obj.w.r = 1-normpdf(obj.f.r,0,shp)/normpdf(0,0,shp);

	% cumulative distribution
	iSr = cumsum(S.r);
	iSr = iSr/iSr(end);
	% make values unique (quick hack)
	iSr = cvec(iSr) + (0:length(iSr)-1)'*1e-12; 
	% quantiles
	f_05 = interp1(iSr,obj.f.r,0.05,'linear');
	f_50 = interp1(iSr,obj.f.r,0.50,'linear');
	f_95 = interp1(iSr,obj.f.r,0.95,'linear');

	% 2D spectral density estimate by periodogram smoothing
	dfr   = frr(1,2);
	nf    = round(sqrt(f_50/dfr));
	S.bar = gaussfilt2(S.clip,nf);

	% smoothing window radius for frequency test
	% TODO no magic numbers
	nf_test  = round(0.25*f_50/dfr);
	nf_test  = max(nf_test,2);

	% restrict test to region containing the upper 20% of spectral energy
	[Ssort,sdx] = sort(S.bar(:));
	iSsort = cumsum(Ssort);
	% TODO no magic numbers
	fdx = (iSsort >= 0.2*iSsort(end));
 	fmsk = false(n);
	fmsk(sdx(fdx)) = true;
%figure(1)
%subplot(2,2,1)
%[Sr_, obj.f.r] = periodogram_radial(S.clip,L);
%plot(obj.f.r,[S.r,Sr_.normalized,obj.w.r*10]);
%vline(fhp)
%plot(obj.f.r,S.r);
%imagesc(fftshift(fmsk))
%imagesc(fftshift(S.clip))

	% exlude spurious low-frequency components from the test
	fmsk = fmsk & (frr > fhp);
%subplot(2,2,2)
%imagesc(fftshift(fmsk))
%pause
	% by symmetry, excluded the symmetric half plane
	fmsk = fmsk & cvec(obj.f.x) >= 0;


	% periodicity test
	%if (numel(bmsk)*obj.opt.ns < obj.opt.n_max_test)
		try
		if (any((bmsk~=bmsk(1,1)),'all'))
			bmsk_ = bmsk;
		else
			bmsk_ = [];
		end
		[p_periodic, stati] = periodogram_test_periodicity_2d(b, ...
								L, nf_test, bmsk_, fmsk, obj.opt.ns);
		fr_periodic = stati.fr_max;
		catch e
			disp(e);
			p_periodic= NaN;
			fr_periodic = NaN;
		end
%	else
%		fprintf('not testing, too large');
%		p_periodic = NaN;
%		fr_periodic = NaN;
%	end

	% determine if pattern is isotropic (spotted, gapped, labyrinthic)
	% or anisotropic (banded)
	% note: presmoothed with Sbar works better than with Shat
	%nt   = pi*n(1);
	%nc   = 2*pi*f_50/dfr;
	nc   = sqrt(f_50/dfr);
	nf_s = ceil(n(1)/nc);
	% more sharper later
	nf_s_ = ceil(1/4*n(1)/nc);
	mode = 'angular'; 
	nf_    = round(2*sqrt((f_50/dfr)));
	Sbar_ = gaussfilt2(S.clip,nf_);
	[isisotropic,stati] = separate_isotropic_from_anisotropic_density(Sbar_,fmsk,L,mode,nf_s);
	angle_deg = stati.angle_deg;
	p_isotropic = stati.p_iso;

	L_eff = effective_mask_size(bmsk,L,-angle_deg);

	% for patterns with known direction, such as computer generated patterns, the direction angle_deg can be specified
	if (isfield('angle',obj.opt) && ~isempty(obj.opt.angle))
		angle_deg = obj.opt.angle;
	%else
		%angle = rad2deg(atan2(fc.yy.bar,fc.xx.bar));
	end

	for field = {'clip','hat','bar'}
		% normalize volume of density to 1
		S.(field{1})   = 2*S.(field{1})/(sum(sum(S.(field{1})))*df(1)*df(2));
		% autocorrelaton
		R.(field{1})   = 0.5*(n(1)*n(2))./(L(1)*L(2))*real(ifft2(S.(field{1})));


		% maximum of the 2d periodogram / density
		[Sc.(field{1}),cdx] = max(S.(field{1}),[],'all');
		%[stat.(field{1}).imax,stat.(field{1}).jmax] = ind2sub(n,cdx);
		[imax,jmax]    = ind2sub(n,cdx);
		% maximum frequency
		fc.rr.(field{1}) = frr(cdx);
		fc.tt.(field{1}) = ftt(cdx);
		fc.xx.(field{1}) = obj.f.x(imax);
		fc.yy.(field{1}) = obj.f.y(jmax);

		% rotate periodogram/density
		S.rot.(field{1}) = fft_rotate(S.(field{1}),-angle_deg);
		R.rot.(field{1}) = fft_rotate(R.(field{1}),-angle_deg);

		% radial periodogram
		[S_,~] = periodogram_radial(S.(field{1}),L);
		S.radial.(field{1}) = S_.normalized;
	
		% angular periodogram
		[S.rot.angular.(field{1}),obj.f.angle] = periodogram_angular(S.rot.(field{1}),L,nf_s_); 


		% density perpendicular to bands
		S.rot.x.(field{1}) = mean(S.rot.(field{1}),2);
		% density parallel to bands
		S.rot.y.(field{1}) = mean(S.rot.(field{1}),1)';

		% normalize the area of density to 1 over the positive half-axis
		S.rot.x.(field{1}) = 2*S.rot.x.(field{1}) / (sum(S.rot.x.(field{1}))*df(1));
		% normalize area under density to 1 over the entire axis
		S.rot.y.(field{1}) = 2*S.rot.y.(field{1}) / (sum(S.rot.y.(field{1}))*df(2));

		% autocorrelation in direction perpendictular to bands
		R.rot.x.(field{1}) = mean(R.rot.(field{1}),2);
		R.rot.x.(field{1}) = R.rot.x.(field{1})/R.rot.x.(field{1})(1);
		%R.rot.x.(field{1}) = 0.5*n(1)/L(1)*real(ifft(S.rot.x.(field{1})));

		% autocorrelation in direction parallel to bands
		%R.rot.y.(field{1}) = 0.5*n(2)/L(1)*real(ifft(S.rot.y.(field{1})));
		R.rot.y.(field{1}) = mean(R.rot.(field{1}),1)';
		R.rot.y.(field{1}) = R.rot.y.(field{1})/R.rot.y.(field{1})(1);

		% radial autocorrelation
		[R.radial.(field{1}),obj.r] = autocorr_radial(R.(field{1}),L);

		% density maxima
		[Sc.radial.(field{1}),id]  = max(S.radial.(field{1}));
		fc.radial.(field{1})       = obj.f.r(id);
		[Sc.x.(field{1}),id]       = max(cvec(S.rot.x.(field{1})).*(obj.f.x>=0));
		fc.x.(field{1})            = obj.f.x(id);
		[Sc.y.(field{1}),id]       = max(cvec(S.rot.y.(field{1})).*(obj.f.y>=0));
		% this should be the first field for anisotropic patterns, but we compute it anyway
		fc.y.(field{1})            = obj.f.y(id);
		[Sc.angular.(field{1}),id] = max(S.rot.angular.(field{1}));
		fc.angular.(field{1})      = obj.f.angle(id);

	end % for field

	fmsk_rot = fft_rotate(fmsk,-angle_deg);

	if (isisotropic)
		regularity = Sc.radial.clip .* fc.radial.clip;
		%[Sc.ref,ldx] = S.radial.clip);
		%lc       = 1./obj.f.r(ldx);
	else
		regularity = Sc.x.clip .* fc.x.clip;
		%[Sc,ldx] = max(Sx);
		%lc       = 1./abs(obj.f.x(ldx));
	end
	

	% store reasults
	stat.L_eff         = L_eff;
	stat.Sc	           = Sc;
	stat.angle_deg     = angle_deg;
	stat.area_msk      = area_msk;
	stat.f_50          = f_50;
	stat.fc            = fc;
	stat.fr_periodic   = fr_periodic;
	stat.isisotropic   = isisotropic;
	stat.fhp           = fhp;
	stat.nf            = nf;
	stat.nf_test       = nf_test;
	stat.p_isotropic   = p_isotropic;
	stat.p_periodic    = p_periodic;
	stat.qfr_05        = f_05;
	stat.qfr_50        = f_50;
	stat.qfr_95        = f_95;
	stat.regularity    = regularity;
	stat.siz           = n;
	stat.centroid      = centroid;

	obj.stat = stat;
	obj.msk.f           = fmsk;
	obj.msk.b_square    = bmsk;
	obj.msk.rot.f       = fmsk_rot;
	obj.R = R;
	obj.S = S;
	obj.f.rr = frr;
	%obj.w = w;
end % analyze_grid

