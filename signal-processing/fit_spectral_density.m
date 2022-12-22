% Tue 20 Jul 22:40:41 CEST 2021
% Karl Kästner, Berlin
%
%% fit spectral densities
%
% function [par,Sp] = fit_spectral_density(fx,S,par,L,form,model,method)
function [par,Sp,cS,stat] = fit_spectral_density(fx,S,par,L,model,method,w,nf,form)
%,isradial)
	if (nargin()<3 || isempty(par))
		error('Initical parameters has to be provided');
		%par = [0.6,1]; 
	end
	if (nargin()<4)
		L = 1;
	end
	if (nargin()<6)
		method = 'ls';
	end
	if (nargin()<7)
		w = 1;
	end
	if (nargin()<8)
		nf = 0;
	end
	if (nargin<9)
		form = 'f';
	end

	S   = cvec(S);
	n   = length(fx);
		
	% restrict always to positive half of the plane
	w = w.*(fx>=0);

	% restrict frequency domain
	wS  = w.*S;

	if (nf > 1)
		wS = ifftshift(trifilt1(fftshift(wS),nf));
	end
	df = 1./L;
	wS = 2*wS/(sum(wS)*df);

	%fdx   = fx>fmin & fx < fmax;
	%IS    = spectral_density_area(fx(fdx),wS(fdx));
	%wS    = wS./IS;
	Sflat = 2*L/n;
	
	switch (method)
	case {'ls','lsrel'}
		opt = optimset();
		opt.Display = 'off';
		%opt.Algorithm = 'levenberg-marquardt';
		lb = sqrt(eps)*zeros(size(par));
		[par,resn,res] = lsqnonlin(@(par) resfun(Sfun(par,nf)),par,lb,[],opt); 
	case {'ll'}
		[par] = fminsearch(@(par) resfun(par),par); 
		resn = 0;
	otherwise
		error('Undefined objective function')
	end

	Sp = Sfun(par,0);

	% compute the r2
	stat.r2 = 1 - wrms(w,Sp-S).^2./wvar(w,S);
	%else
	%	stat.r2 = 1 - wrms(abs(fx(fdx)),wSp(fdx)-wS(fdx)).^2./wvar(abs(fx(fdx)),wS(fdx));
	%end
	stat.resn = resn;
	

	% unfiltered
	%Sp = Sfun(par,0);
	%wwSp = ifftshift(trifilt1(fftshift(wSp),nf));
	% bias corrected density estimate
	%if (0)
	%cS = wSp.*Sp./wwSp;
	%else
	cS = Sp;
	%end

	function Sp = Sfun(par,nf)
		switch (model)
		case {'bandpass-continuous'}
			Sp = spectral_density_bandpass_continuous(fx,par(1),par(2));
		case {'bandpass-continuous-mean'}
			Sbp = spectral_density_bandpass_continuous(fx,par(1),par(2));
			IS_    = spectral_density_area(fx,S);
			S_  = S/IS_;
			p0  = 2*(S_(2)-Sbp(2))-(S_(3)-Sbp(3));
			p0  = max(p0,0);
			Sm  = exp(-par(3)*abs(fx));
			Sp  = p0*Sm + (1-p0/par(3))*Sbp;                                                  
		case {'bandpass-continuous-mean-old'}
			Sp_bp   = spectral_density_bandpass_continuous(fx,par(1),par(2));
			Sp_mean = exp(-par(3)*abs(fx));
			Sp      = par(4)*Sp_bp + (1-par(4))*Sp_mean;
		case {'bandpass-discrete'}
			dx = L/n;
			Sp = spectral_density_bandpass_discrete(fx,par(1),par(2),dx,form);
		case {'lorentzian'}
			Sp = spectral_density_lorentzian(fx,par(1),par(2));
		case {'brownian-phase'}
			Sp = spectral_density_brownian_phase(fx,par(1),par(2));
		case {'brownian-phase-mean'}
			Sbm = spectral_density_brownian_phase(fx,par(1),par(2));
			IS_    = spectral_density_area(fx,S);
			S_  = S/IS_;
			p0  = 2*(S_(2)-Sbm(2))-(S_(3)-Sbm(3));
			p0  = max(p0,0);
			Sm  = exp(-par(3)*abs(fx));
			Sp  = p0*Sm + (1-p0/par(3))*Sbm;                                                   
		case {'brownian-phase-mean-old'}
			Spbm   = spectral_density_brownian_phase(fx,par(1),par(2));
			Spbm   = Spbm/sum(Spbm);
			Spmean = spectral_density_brownian_phase_across(fx,par(3));
			Spmean =Spmean/sum(Spmean);
			Sp = par(4)*Spbm + (1-par(4))*Spmean;
		case {'brownian-phase-across'}
			Sp = spectral_density_brownian_phase_across(fx,par(1));
		case {'lognormal'}
			Sp = lognpdf(abs(fx),par(1),par(2));
			Sp = lognpdf(abs(fx),par(1),par(2));
		otherwise
			disp(model);
			error('Undefined spectral density function')
		end
		Sp = real(Sp);
		%if (nf > 1)
		%	Sp = ifftshift(trifilt1(fftshift(Sp),nf));
		%end
		ISp = spectral_density_area(fx,Sp);
		Sp = Sp/ISp;
	end

	function res = resfun(Sp)
		% restrict frequency domain (to exclude for example spurious low frequency lobes)
		wSp = w.*Sp;

		% smooth, this is necessary, to assure continuity for very spiky densities of nearly periodic patterns
		if (nf > 1)
			wSp = ifftshift(trifilt1(fftshift(wSp),nf));
		end
		% normalize
		wSp = 2*wSp/(sum(wSp)*df);

		switch (method)
		case {'lsrel'}
			res = (wS-wSp)./(wSp+wSflat);
			%res   = (wwS-wwSp)./(wSp(fdx)+Sflat);
			%res   = (wS(fdx)-wSp(fdx))./(wSp(fdx)+Sflat);
		case {'ls'}
			res = wSp - wS;
			%if (isradial)
			%	res = abs(fx(fdx))
			%end
		case {'ll'}
			res = log(wS./(wSp+Sflat));
			%res    = (log(Sp+Sflat) + S./(Sp+Sflat));
			% note, this assumes that fx(1) = 0
			res(1) = 0;
			res    = sum(res(fdx));
		end
	end
end % fit_spectral_density

