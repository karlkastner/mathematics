% Tue 20 Jul 22:40:41 CEST 2021
% Karl Kästner, Berlin
%
%% fit (spectral) densities
%
% output :	par   : parameters of fitted density
%		Sp    : fitted density
%		stat  : statistics of fist
% input  :	f     : frequencies, at which the input density/periodogram is sampled
%		       only the positive half-axis (f>=0) is used for the fit
%		S     : 1D density or periodogram
%		par0  : initial parameter for fit
%		density_model : density model
%		method: objective to optimize
%		w     : weight certain frequency components,
%		       can be used to suppress frequency components when w(f) = 0
%		nf    : radius of smoothing window (for least-squares fits)
%
% method:
%	mise-cramer	: minimize ( cdf(wS) - cdf(whatS))
%			recommended default
%			- robust against added noise
%			- converges even when S is nearly discrete
%			- does not require smoothing
%
%	log-likelyhood : maximize log(whatS/wS(p))
%			 not recommended, as not robust against added noise to hat S
%			 can fail to converge when distribution becomes discrete
%
%	least square : minimize least sqares log(bar (w hat S - S))
%		         not recommended
%			 - reqieres smoothing
%			 - robust against added noise
%			 - can fail to converge when distribution becomes discrete
%
% function [par,Sp] = fit_spectral_density(f,S,w,density_model,par0,method,nf)
function [par,Sp,stat] = fit_spectral_density(f,S,w,density_model,par0,method,nf)
	if (nargin()<3 || isempty(par0))
		error('Initical parameters has to be provided');
		%par = [0.6,1]; 
	end
	if (nargin()<4)
		L = 1;
	end
	if (nargin()<6)
		method = 'mise-cramer';
	end
	if (nargin()<7)
		w = 1;
	end
	if (nargin()<8)
		nf = 3;
	end

	S   = cvec(S);
	n   = length(f);
	
	% store initial f	
	f_ = f;
	% restrict fit to positive half of the plane
	fdx = f>0;
	f = f(fdx);
	S = S(fdx);
	if (~isscalar(w))
		w = w(fdx);
	end
	IS_    = spectral_density_area(f,S);
	S_  = S/IS_;
	
%	w = w.*(f>=0);

	% weigh input distribution
	% (mainly for suppression of spurious-low frequency components)
	wS  = w.*S;

	%df = 1./L;
	% renormalize
	wS = wS/sum(mid(wS).*diff(f));

	%fdx   = f>fmin & f < fmax;
	%IS    = spectral_density_area(f(fdx),wS(fdx));
	%wS    = wS./IS;
	%Sflat = 2*L/n;
	
	switch (method)
	case {'mise-cramer','mc'}
		iwS = cumint_trapezoidal(f,wS);
		opt = optimset();
		opt.Display = 'off';
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) res_mise_cramer(Sfun(par)), par0,[],[],opt); 
	case {'least-sqares','ls'}
		opt = optimset();
		opt.Display = 'off';
		%opt.Algorithm = 'levenberg-marquardt';
		lb = sqrt(eps)*zeros(size(par0));
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) res_least_squares(Sfun(par)),par0,lb,[],opt); 
	case {'log-likelihood','ll'}
		par  = fminsearch(@(par) res_log_likelyhood(par),par0); 
		resn = 0;
	otherwise
		error('Undefined objective function')
	end

	% compute the r2 of windowed smoothed density
	% (n.b. : this is not a good measure, the mise-cramer distance is superior)
	Sp = Sfun(par);
	stat.r2 = 1 - wrms(w,Sp-S).^2./wvar(w,S);

	% predict density without smoothing
	nf = 0;
	f = f_;
	Sp = Sfun(par);

	% unfiltered
	%Sp = Sfun(par,0);
	%wwSp = ifftshift(trifilt1(fftshift(wSp),nf));
	% bias corrected density estimate
	%if (0)
	%cS = wSp.*Sp./wwSp;
	%else
	%cS = Sp;
	%end

	function Sp = Sfun(par)
		switch (density_model)
		case {'bandpass-continuous'}
			Sp = spectral_density_bandpass_continuous(f,par(1),par(2));
		case {'bandpass-continuous-mean'}
			Sbp = spectral_density_bandpass_continuous(f,par(1),par(2));
			p0  = 2*(S_(2)-Sbp(2))-(S_(3)-Sbp(3));
			p0  = max(p0,0);
			Sm  = exp(-par(3)*abs(f));
			Sp  = p0*Sm + (1-p0/par(3))*Sbp;                                                  
		case {'bandpass-continuous-mean-old'}
			Sp_bp   = spectral_density_bandpass_continuous(f,par(1),par(2));
			Sp_mean = exp(-par(3)*abs(f));
			Sp      = par(4)*Sp_bp + (1-par(4))*Sp_mean;
		case {'bandpass-discrete'}
			dx = L/n;
			Sp = spectral_density_bandpass_discrete(f,par(1),par(2),dx);
		case {'lorentzian'}
			Sp = spectral_density_lorentzian(f,par(1),par(2));
		case {'brownian-phase'}
			Sp = spectral_density_brownian_phase(f,par(1),par(2));
		case {'brownian-phase-mean'}
			Sbm = spectral_density_brownian_phase(f,par(1),par(2));
			%IS_    = spectral_density_area(f,S);
			%S_  = S/IS_;
			p0  = 2*(S_(2)-Sbm(2))-(S_(3)-Sbm(3));
			p0  = max(p0,0);
			Sm  = exp(-par(3)*abs(f));
			Sp  = p0*Sm + (1-p0/par(3))*Sbm;                                                   
		case {'brownian-phase-mean-old'}
			Spbm   = spectral_density_brownian_phase(f,par(1),par(2));
			Spbm   = Spbm/sum(Spbm);
			Spmean = spectral_density_brownian_phase_across(f,par(3));
			Spmean =Spmean/sum(Spmean);
			Sp = par(4)*Spbm + (1-par(4))*Spmean;
		case {'brownian-phase-across'}
			Sp = spectral_density_brownian_phase_across(f,par(1));
		case {'lognormal'}
			Sp = lognpdf(abs(f),par(1),par(2));
			Sp = lognpdf(abs(f),par(1),par(2));
		otherwise
			disp(density_model);
			error('Undefined spectral density model')
		end
		Sp = real(Sp);
		%if (nf > 1)
		%	Sp = ifftshift(trifilt1(fftshift(Sp),nf));
		%end
		ISp = spectral_density_area(f,Sp);
		Sp = Sp/ISp;
	end

	% note that the residual does not need to squared and summed
	% lsqnonlin does this implicitely and more efficient when doing so
	function res = res_least_squares(Sp)
		% weigh
		wSp = w.*Sp;
		% renormalize
		wSp = wSp/sum(mid(wSp).*diff(f));	
		% residual
		res = wSp - wS;
		% smooth
		res = trifilt1(res,nf);	
	end

	function res = res_log_likelyhood(Sp)
		% weigh
		wSp = w.*Sp;
		% renormalize
		wSp = wSp/sum(mid(wSp).*diff(f));	
		% log-likelyhood
		res   = log(wS./wSp);
		%res   = (log(Sp+Sflat) + S./(Sp+Sflat));
		res(wSp==0) = 0;
		res    = sum(res);
	end
	
	function res = res_mise_cramer(Sp)
		% weigh
		wSp = w.*Sp;
		% renormalize
		wSp = wSp/sum(mid(wSp).*diff(f));
		% cdf of weighed distribution
		iwSp = cumint_trapezoidal(f,wSp);
		res = iwSp - iwS;
	end
end % fit_spectral_density

