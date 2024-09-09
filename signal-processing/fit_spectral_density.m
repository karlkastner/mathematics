% Tue 20 Jul 22:40:41 CEST 2021
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% fit spectral densities (probability distributions)
%
% output :	par   : parameters of fitted density
%		Sp    : fitted density
%		stat  : statistics of fist
% input  :	f     : frequencies, at which the input density/periodogram is sampled
%		       only the positive half-axis (f>=0) is used for the fit
%		S     : 1D density or periodogram
%		par0  : initial parameter for fit
%		Sfun  : density model function (propobility distribution)
%		distance: objective to optimize
%		w     : weight certain frequency components,
%		       can be used to suppress frequency components when w(f) = 0
%		nf    : radius of smoothing window (for least-squares fits)
%
% distance: objective to be minimized
%
%	mise-cramer	: minimize ( cdf(wS) - cdf(whatS))
%			recommended default
%			- robust against added noise
%			- converges even when S is nearly discrete
%			- does not require smoothing
%
%	log-likelihood : maximize log(whatS/wS(p))
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
function [par,Sp,stat] = fit_spectral_density(f,S,w,Sfun,par0,distance,nf,lb)
	if (nargin()<3)
		w = 1;
	end
	if (nargin()<4 || isempty(par0))
		error('Initial parameters has to be provided');
	end
	if (nargin()<6)
		distance = 'hellinger';
	end
	if (nargin()<7||isempty(nf))
		nf = 1;
	end
	if (nargin()<8)
		lb = [];
	end
	% for lsqnonlin
	S = double(S);

	% exclude negative part, as it is redundant
	% note that ther is otherwise a jump without shifting the fft
	fdx = f>=0;
	f_   = cvec(f(fdx));
	S_   = cvec(S(fdx));
	df_  = diff(f_);
	sqrt_df = sqrt(df_);

	%n   = length(f);

	% store initial f
	%f = f;

	if (isa(w,'function_handle'))
		w = w(f);
		w  = cvec(w);
	else
		if (~isscalar(w))
			w  = cvec(w);
			w_ = w(fdx);
		end
	end
	w  = double(w);
	w_ = double(w_);

	% restrict fit to positive half of the plane
	%fdx  = f>0;
	%fpos = f(fdx);
	%Spos = S(fdx);
	%if (~isscalar(w))
	%	w = w(fdx);
	%end

	% normalize
	%int_S = spectral_density_area(f,S);
	%S     = S/int_S;

	% weigh input distribution
	% (mainly for suppression of spurious-low frequency components)
	wS_  = w_.*S_;

	% normalize
	int_wS_ = spectral_density_area(f_,wS_);
	wS_ = wS_/int_wS_;
	%sum(mid(wS).*diff(f));

	iwS_ = cumint_trapezoidal(f_,wS_);

	par0 = double(par0);
	wSp0 = wSfun(par0);
	wSp0 = double(wSp0);

%	wSfun_ = @(x) double(wSfun(x));

	switch (distance)
	case {'cdf_l1'}
		opt = optimset();
		opt.Display = 'off';
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) res_cdf_l1(wSfun(par)), par0,lb,[],opt);
	case {'cdf_l2','cramer-von-mises','mvc'}
		opt = optimset();
		opt.Display = 'off';
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) res_mise_cramer(wSfun(par)), par0,lb,[],opt);
	case {'absolute-difference'} % note that this is equal to minimizing the maximum distance of the cdf

	case {'least-squares','ls'}
		opt = optimset();
		opt.Display = 'off';
		%opt.Algorithm = 'levenberg-marquardt';
		%lb = sqrt(eps)*zeros(size(par0));
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) res_least_squares(wSfun(par)),par0,lb,[],opt);
	case {'hellinger'}
		sqrt_wS = sqrt(wS_);
		% TODO normalize by L
		opt.Display = 'off';
		[par,stat.resn,res,stat.exitflag] = lsqnonlin(@(par) sqrt(wSfun(par)) - sqrt_wS,par0,lb,[],opt);
	case {'log-likelihood','ll'}
		par  = fminsearch(@(par) log_likelihood(par),par0);
		resn = 0;
	otherwise
		error('Undefined objective function')
	end

	parc   = num2cell(par);
	Sp     = Sfun(f,parc{:});
	int_Sp = spectral_density_area(f,Sp);
	Sp     = Sp/int_Sp;
	Sp_    = Sfun(f_,parc{:});
	int_Sp_ = spectral_density_area(f,Sp_);
	Sp_     = Sp_/int_Sp_;
	wSp_   = wSfun(par);

	% goodness of fit measures
	% root mean square error
	res_sqrt_df          = res_least_squares(wSp_);
	% rmse = sqrt(1/(f_max-0) int_0^{f_max} res.^2 df)
	%stat.goodness.rmse   = sqrt(1./max(f)*sum(res_sqrt_df.*res_sqrt_df));
	%stat.goodness.r2     = 1.0 - stat.goodness.rmse.^2./wvar(w_.*inner2outer(df_),S_);
	stat.goodness.hd      = hellinger_distance(S_,Sp_,df_(1),w_);
	stat.goodness.r2      = 1 - 2*stat.goodness.hd;

	% l2 in difference of cdf von Mise - Cramer distance
	stat.goodness.mise_cramer = res_mise_cramer(wSp_);

	% l1 of difference in cdf
	stat.goodness.cdf_l1 = res_cdf_l1(wSp_);

	% log-likelihood
	stat.goodness.log_likelihood = log_likelihood(wSp_);

	% compute the r2 of windowed smoothed density
	% (n.b. : this is not a good measure, the mise-cramer distance is superior)
	%Sp = Sfun(par);
	%stat.r2 = 1 - wrms(w,Sp-S).^2./wvar(w,S);

%	function weight(par)
%		parc = num2cell(par);
%		Sp_  = Sfun(f_,parc{:});
%	end

	function wSp_ = wSfun(par)
		parc = num2cell(par);
		Sp_  = Sfun(f_,parc{:});
		% weigh
		wSp_ = w_.*Sp_;
		% normalize
		int_wSp_ = spectral_density_area(f_,wSp_);
		wSp_ = wSp_/int_wSp_;
	end

	% note that the residual does not need to squared and summed
	% lsqnonlin does this implicitely and more efficient when doing so
	function res = res_least_squares(wSp_)
		% residual
		res = wSp_ - wS_;
		res = mid(res).*sqrt_df;
	end

	function ll = log_likelihood(wSp_)
		% weigh
		%wSp = w.*Sp;
		% renormalize
		%wSp = wSp/sum(mid(wSp).*diff(f));
		% log-likelihood
		ll   = log(wS_./wSp_);
		%res   = (log(Sp+Sflat) + S./(Sp+Sflat));
		ll(wSp_==0) = 0;
		ll  = sum(mid(ll).*df_);
	end

	function res = res_mise_cramer(wSp_)
		% weigh
		%wSp = w.*Sp;
		% renormalize
		%wSp = wSp/sum(mid(wSp).*diff(f));
		% cdf of weighed distribution
		iwSp_ = cumint_trapezoidal(f_,wSp_);
		res  = iwSp_ - iwS_;
%		res  = mid(res).*sqrt_df;
%		res  = mid(res); %.*sqrt_df;
%		res = sum(res.^2);
		res = sum(mid(res.*res).*df_);
	end

	function res = res_cdf_l1(wSp_)
		% weigh
		%wSp = w.*Sp;
		% renormalize
		%wSp = wSp/sum(mid(wSp).*diff(f));
		% cdf of weighed distribution
		iwSp_ = cumint_trapezoidal(f_,wSp_);
		res  = iwSp_ - iwS_;
%		res  = mid(res).*sqrt_df;
%		res  = mid(res); %.*sqrt_df;
%		res = sum(res.^2);
		res = sum(mid(abs(res)).*df_);
	end
end % fit_spectral_density

