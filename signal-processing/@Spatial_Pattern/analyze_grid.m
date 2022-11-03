% Wed 18 May 13:50:47 CEST 2022
function analyze_grid(obj)
	% output
	S    = struct();
	R    = struct();
	stat = struct();

	n = obj.n;
	L = obj.L;
		
	b   = obj.b;
	bmsk = obj.msk;
	if (isempty(bmsk))
		bmsk = true(size(b));
	end

	% images might be of integer type and have to be converted 
	b    = double(b);
	bmsk = bmsk > 0;

	% fraction of ground coverged
	% there is a matlab bug that double images need to be scaled
	bscaled = (b - min(b(bmsk)))/(max(b(bmsk))-min(b(bmsk)));

	thresh_b      = graythresh(bscaled(bmsk));
	stat.coverage = sum(bscaled(bmsk)<thresh_b)/sum(bmsk(:));
	
	% resample, to make dx identical to dy
	% (depending on the grid projection, this might not be the case)
	dxy = L./n;
	n_ = round(L./min(dxy));
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
	df   = 1./L;
	n    = [n,n];
	b(end+1:nmax,:)   = 0;
	b(:,end+1:nmax)   = 0;
	bmsk(end+1:nmax,:) = false;
	bmsk(:,end+1:nmax) = false;
	n = size(b);

	% grid in real space
	% TODO x an y do not match original figure when the figure was resampled or or resized
	obj.x   = linspace(-L(1)/2,L(1)/2,n(1))';
	obj.y   = linspace(-L(2)/2,L(2)/2,n(2))';
	rr      = hypot(obj.x',obj.y);

	% grid in frequency space
	[obj.f.x,obj.f.y,obj.f.rr,f.tt] = fourier_axis_2d(L,n);

	% TODO smooth mask with border

	% remove offset from 0
	wbmsk = double(bmsk);
	%b_ = (b - mean(flat(wbmsk.*b)));
	b_ = (b -  wmean(wbmsk(:),b(:)));

	% normalize
	b_ = b_/wrms(wbmsk(:),b_(:));

	% periodogram
	% TODO use function
	% TODO padd zeros to reach similar L for analysis for all patterns
	% TODO iterate with clipping?
	S.raw  = abs(fft2(wbmsk.*b_)).^2;

	% mean frequency component
	fr_mean = wmean(S.raw(:),obj.f.rr(:))


	f_min = obj.opt.fminscale*fr_mean;
	f_max = obj.opt.fmaxscale*fr_mean;

	% remove spurious low-frequency components by lowpass filtering (truncation)
%	S.raw = fftshift(S.raw);

	% windows for frequency range of interest
	w.frr = ones(size(obj.f.rr));
	w.f.x  = ones(size(obj.f.x));
	w.f.y  = ones(size(obj.f.y));

	% mask spurious low-frequency components
	s = sqrt(-log(0.25));
	if (~isempty(f_min) && f_min > 0)
		w.frr = w.frr .* (1-normpdf(s*obj.f.rr/f_min)/normpdf(0));
		w.f.x  = w.f.x  .* (1-normpdf(s*obj.f.x/f_min)/normpdf(0));
		w.f.y  = w.f.y  .* (1-normpdf(s*obj.f.y/f_min)/normpdf(0));
	end
	% mask high pass frequencies beyond range of interest
	if (~isempty(f_max) && f_max < inf)
		w.frr = w.frr .* ( normpdf(s*obj.f.rr/f_max)/normpdf(0));
		w.f.x  = w.f.x  .* ( normpdf(s*obj.f.x/f_max)/normpdf(0));
		w.f.y  = w.f.y  .* ( normpdf(s*obj.f.y/f_max)/normpdf(0));
	end
	rmax        = 5/fr_mean;
	w.fsmooth   = radial_window(ifftshift(rr),rmax);

	% restrict to frequency range of interest
	S.clip = w.frr.*S.raw;

	% spectral density through smoothing
	S.gauss = real(fft2(w.fsmooth.*ifft2(S.clip)));

	for field = {'raw','clip','gauss'}
		% normalize volume of density to 1
		S.(field{1})   = 2*S.(field{1})/(sum(S.(field{1})(:))*df(1)*df(2));
		% autocorrelaton
		R.(field{1})   = 0.5*real(ifft2(S.(field{1})));
		% normalize autocorrelation to 1
		R.(field{1}) = R.(field{1})/R.(field{1})(1,1);
		% radial periodogram
		[Smu, S.radial.(field{1}),kSr,obj.f.r] = periodogram_radial(S.(field{1}),L);
		% radial autocorrelation
		[R.radial.(field{1}),obj.r] = autocorr_radial(R.(field{1}),L);
		% maximum of the periodogram / density
		[Sc.(field{1}),cdx] = max(S.(field{1}),[],'all');
		%[stat.(field{1}).imax,stat.(field{1}).jmax] = ind2sub(n,cdx);
		[imax,jmax]    = ind2sub(n,cdx);
		% maximum frequency
		fc.rr.(field{1}) = obj.f.rr(cdx);
		fc.tt.(field{1}) = f.tt(cdx);
		fc.xx.(field{1}) = obj.f.x(imax);
		fc.yy.(field{1}) = obj.f.y(jmax);
	end % for field

	w.f.r  = ones(size(obj.f.r));
	w.f.r  = w.f.r  .* (1 - normpdf(s*obj.f.r/f_min)/normpdf(0));
	w.f.r  = w.f.r  .* (normpdf(s*obj.f.r/f_max)/normpdf(0));

	%stat.nf = 10;
	% TODO dynamic
	% TODO, this is related to the gaussiam density estimate (also smoothing)
	stat.nf_test = 10;
	% = max(7,ceil(0.3*fcr./df));

	% periodicity test
	% TODO the periodicity test has a deviating threshold when the bmask is not everywhere 1
	[stat.pt, stat_, ratio, stat.qq] = periodogram_test_periodicity_2d(b, obj.L, stat.nf_test, bmsk, f_min, f_max);

	%stat.f_radial = f_radial;

	% this is computed anyway for all patterns, but only meaningfull for banded patterns
	if (1)
		% for patterns with known direction, such as computer generated patterns, the direction angle can be specified
		if (isfield('angle',obj.opt) && ~isempty(obj.opt.angle))
			angle = obj.opt.angle;
		else
			angle = rad2deg(atan2(fc.yy.gauss,fc.xx.gauss));
		end

		for field = {'raw','clip','gauss'}
			% rotate periodogram/density
			S.rot.(field{1})   = ifftshift(imrotate(fftshift(S.(field{1})),  -angle,'bilinear','crop'));
			%S.rot.gauss = ifftshift(imrotate(fftshift(S.gauss),-stat.angle,'bilinear','crop'));
			%R.rot.(field{1})   = ifftshift(imrotate(fftshift(R.(field{1})),  -angle,'bilinear','crop'));

			% density perpendicular to bands
			S.rot.x.(field{1}) = mean(S.rot.(field{1}),2);
			% density parallel to bands
			S.rot.y.(field{1}) = mean(S.rot.(field{1}),1)';

			% normalize area of density to 1
			S.rot.x.(field{1}) = 2*S.rot.x.(field{1}) / (sum(S.rot.x.(field{1}))*df(1));
			S.rot.y.(field{1}) = 2*S.rot.y.(field{1}) / (sum(S.rot.y.(field{1}))*df(2));

			% autocorrelation in direction perpendictular to bands
			R.rot.x.(field{1}) = 2*real(ifft(S.rot.x.(field{1})));
			% autocorrelation in direction parallel to bands
			R.rot.y.(field{1}) = 2*real(ifft(S.rot.y.(field{1})));

			% density maxima
			%fdx   = (obj.f.r >= f_min) & (obj.fr<= f_max);
			[Sc.r.(field{1}),id]   = max(S.radial.(field{1}));
			fc.r.(field{1})   = obj.f.r(id);
			%fdx   = (obj.f.x>= f_min) & (obj.f.x<= f_max);
			[Sc.x.(field{1}),id]   = max(cvec(S.rot.x.(field{1})).*(obj.f.x>=0));
			fc.x.(field{1})   = obj.f.x(id);
			%fdx   = (obj.f.y>= f_min) & (obj.f.y<= f_max);
			[Sc.y.(field{1}),id]   = max(cvec(S.rot.y.(field{1})).*(obj.f.y>=0));
			fc.y.(field{1})   = obj.f.y(id);
		end

		% fit parametric density models in direction perpendicular to bands
		try

			% without mean component
			[par(1),par(2)] = spectral_density_brownian_phase_mode2par(fc.rr.gauss,Sc.x.gauss);
			par(2) = 0.1;
			nf  = 3;
			[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.raw,par,L(1),'brownian-phase','ls',w.f.x,nf);
			Sfit = Sfit/sum(Sfit)*sum(S.rot.x.raw);
			S.rot.x.brownian_phase = Sfit;
			stat.fit.x.brownian_phase.par  = par;
			stat.fit.x.brownian_phase.stat = fitstat;

			par(3) = 5./fc.rr.gauss;
			[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.raw,par,L(1),'brownian-phase-mean','ls',w.f.x,nf);
			Sfit = Sfit/sum(Sfit)*sum(S.rot.x.raw);
			S.rot.x.brownian_phase_mean = Sfit;
			stat.fit.x.brownian_phase_mean.par  = par;
			stat.fit.x.brownian_phase_mean.stat = fitstat;

			% without mean component
			par = [fc.rr.gauss,20];
			[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.raw,par,L(1),'bandpass-continuous','ls',w.f.x,nf);
			Sfit = Sfit/sum(Sfit)*sum(S.rot.x.raw);
			S.rot.x.bandpass         = Sfit;
			stat.fit.x.bandpass.par  = par;
			stat.fit.x.bandpass.stat = fitstat;

			% with mean
			par(3) = [10./fc.rr.gauss];
			[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.raw,par,L(1),'bandpass-continuous-mean','ls',w.f.x,nf);
			Sfit = Sfit/sum(Sfit)*sum(S.rot.x.raw);
			S.rot.x.bandpass_mean = Sfit;
			stat.fit.x.bandpass_mean.par  = par;
			stat.fit.x.bandpass_mean.stat = fitstat;

			%cS_ = cS_/sum(cS_)*sum(cS);
			%plot(f.rc/fc(kdx),cS_/sqrt(cSc)*fc(kdx),'linewidth',1.5);
%		catch e
%			e
%		end

		% fit parametric density models in the direction perpendicular to the bands
%		try
			par3 = 1;
			fa = obj.f.y;
			[par3,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.y,S.rot.y.raw,par3,L(2),'brownian-phase-across','ls',w.f.y);
			S.rot.y.brownian_phase_across = Sfit;
			stat.fit.y.brownian_phase_across.par  = par3;
			stat.fit.y.brownian_phase_across.stat = fitstat;
			%aS_ = aS_/(sum(aS_)*(fa(2)-fa(1)));
%		catch e
%			e
%		end
	
	% fit parametric density models to the radial density
		nf = 1;
		par = [fc.rr.gauss,0.1];
		Lr = 1./(obj.f.r(2)-obj.f.r(1));
		[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.raw,par,Lr,'brownian-phase','ls',w.f.r,nf);
		%cS_ = cS_/sum(cS_)*sum(cS);
		S.radial.brownian_phase = Sfit;
		stat.fit.radial.brownian_phase.par  = par;
		stat.fit.radial.brownian_phase.stat = fitstat;

%	try
		par(3) = 5./fc.rr.gauss;
		[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.raw,par,Lr,'brownian-phase-mean','ls',w.f.r,nf);
		%cS_ = cS_/sum(cS_)*sum(cS);
		S.radial.brownian_phase_mean = Sfit;
		stat.fit.radial.brownian_phase_mean.par  = par;
		stat.fit.radial.brownian_phase_mean.stat = fitstat;

		par = [fc.rr.gauss,10];
		nf = 1;
		[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.raw,par,Lr,'bandpass-continuous','ls',w.f.r,nf);
		S.radial.bandpass = Sfit;
		stat.fit.radial.bandpass.par  = par;
		stat.fit.radial.bandpass.stat = fitstat;
%	try
		par(3) = 10./fc.rr.gauss;
		[par,Sfit_,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.raw,par,Lr,'bandpass-continuous-mean','ls',w.f.r,nf);
		S.radial.bandpass_mean = Sfit;
		stat.fit.radial.bandpass_mean.par  = par;
		stat.fit.radial.bandpass_mean.stat = fitstat;
		%plot(f.rc/fc(kdx),cS_/sqrt(cSc)*fc(kdx),'linewidth',1.5);
%	catch e
%		e
%	end
	% central frequency, radius and angle
	end

	stat.Sc	     = Sc;
	stat.fc      = fc;
	stat.frr_mean = fr_mean;
	%stat.mask     = mask;

	obj.R = R;
	obj.S = S;
	obj.w = w;
	obj.stat = stat;
end % analyze_grid

