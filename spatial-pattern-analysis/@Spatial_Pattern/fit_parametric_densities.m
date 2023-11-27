% Wed 18 May 13:50:47 CEST 2022
% Karl Kastner, Berlin
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
%% fit parametric spectral densities to the empirical density
function fit_parametric_densities(obj)
	S    = obj.S;
	stat = obj.stat;
	fc = obj.stat.fc;
	Sc = obj.stat.Sc;

	% TODO from obj
	nf = [];

	%
	% fit parametric density models in direction perpendicular to bands
	% this is computed anyway for all patterns, but only meaningfull for banded patterns
	%

	% phase drift
	[par0(1),par0(2)] = phase_drift_pdf_mode2par(fc.x.bar,Sc.x.bar);
	Sfun = @phase_drift_pdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,Sfun,par0,obj.opt.objective,nf,[0,0]);
	S.rot.x.phase_drift = Sfit;
	stat.fit.x.phase_drift.par  = par;
	stat.fit.x.phase_drift.stat = fitstat;

	% bandpass
	par0 = [fc.rr.bar, bandpass1d_continuous_pdf_max2par(fc.x.bar,Sc.x.bar,10)];
	Sfun = @bandpass1d_continuous_pdf;	
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,Sfun,par0,obj.opt.objective,nf,[0,0.5]);
	S.rot.x.bandpass         = Sfit;
	stat.fit.x.bandpass.par  = par;
	stat.fit.x.bandpass.stat = fitstat;

	% log-normal
	[par0(1),par0(2)] = logn_mode2par(fc.x.bar,Sc.x.bar);
	[par0(1),par0(2)] = logn_mode2par(fc.x.hat,Sc.x.hat);
	Sfun = @lognpdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,Sfun,par0,obj.opt.objective,nf,[-inf,0]);
	S.rot.x.logn  = Sfit;
	stat.fit.x.logn.par  = par;
	stat.fit.x.logn.stat = fitstat;

	% gamma
	[par0(1),par0(2)] = gamma_mode2par(fc.x.hp,Sc.x.hp);
	par0(1) = max(par0(1),1);
	Sfun = @gampdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,Sfun,par0,obj.opt.objective,nf);
	S.rot.x.gamma  = Sfit;
	stat.fit.x.gamma.par  = par;
	stat.fit.x.gamma.stat = fitstat;

	% TODO wrapped-normal

	% fit parametric density models in the direction perpendicular to the bands
	fy = obj.f.y;
	fdx = (fy>=0);
	par0 = phase_drift_parallel_pdf_mode2par(S.rot.y.hat(1));
	Sfun = @phase_drift_parallel_pdf;
	[pary,Sfit,fitstat] = fit_spectral_density(obj.f.y,S.rot.y.hat,obj.w.y,Sfun,par0,obj.opt.objective,nf);
	S.rot.y.phase_drift_parallel = Sfit;
	stat.fit.y.phase_drift_parallel.par  = pary;
	stat.fit.y.phase_drift_parallel.stat = fitstat;


	% flat density of white noise
	% TODO weight
	S.fit.x.white = mean(S.rot.x.hat)*ones(size(S.rot.x.hat));
	stat.fit.x.white.stat.goodness.r2 = 1 - rms(S.rot.x.hat - S.fit.x.white).^2./var(S.rot.x.hat);

	% periodic
	flag = obj.f.x>=0;
	[mv,mdx] = max(S.rot.x.hat.*flag); 
	S.fit.x.periodic = zeros(size(S.rot.x.hat));
	S.fit.x.periodic(mdx) = sum(S.rot.x.hat.*flag);
	stat.fit.x.periodic.stat.goodness.r2 = 1 - rms(S.rot.x.hat(flag) - S.fit.x.periodic(flag)).^2./var(S.rot.x.hat(flag));


	% fit parametric density models to the radial density

	% phase drift
	par0 = [];
	[par0(1),par0(2)] = phase_drift_pdf_mode2par(fc.radial.bar,Sc.radial.bar);
	Sfun = @phase_drift_pdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,Sfun,par0,obj.opt.objective,nf);
	S.radial.phase_drift = Sfit;
	stat.fit.radial.phase_drift.par  = par;
	stat.fit.radial.phase_drift.stat = fitstat;

	% bandpass
	par0 = [fc.radial.bar, bandpass1d_continuous_pdf_max2par(fc.radial.bar,Sc.radial.bar,10)];
	Sfun = @bandpass1d_continuous_pdf;	
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,Sfun,par0,obj.opt.objective,nf,[0,0.5]);
	S.radial.bandpass = Sfit;
	stat.fit.radial.bandpass.par  = par;
	stat.fit.radial.bandpass.stat = fitstat;

	% log-normal
	[par0(1),par0(2)] = logn_mode2par(fc.radial.hat,Sc.radial.hat);
	Sfun = @lognpdf;

	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,Sfun,par0,obj.opt.objective,nf,[-inf,0]);
	S.radial.logn  = Sfit;
	stat.fit.radial.logn.par  = par;
	stat.fit.radial.logn.stat = fitstat;

	% gamma
	[par0(1),par0(2)] = gamma_mode2par(fc.radial.bar,Sc.radial.bar);
	Sfun = @gampdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,Sfun,par0,obj.opt.objective,nf);
	S.radial.gamma  = Sfit;
	stat.fit.radial.gamma.par  = par;
	stat.fit.radial.gamma.stat = fitstat;

	% flat density of white noise
	% TODO weight
	S.fit.radial.white = mean(S.radial.hat)*ones(size(S.radial.hat));
	stat.fit.radial.white.stat.goodness.r2 = 1 - rms(S.radial.hat - S.fit.radial.white).^2./var(S.radial.hat)

	% periodic
	flag = obj.f.r>=0;
	[mv,mdx] = max(S.radial.hat.*flag); 
	S.fit.radial.periodic = zeros(size(S.radial.hat));
	S.fit.radial.periodic(mdx) = sum(S.radial.hat.*flag);
	stat.fit.radial.periodic.stat.goodness.r2 = 1 - rms(S.radial.hat(flag) - S.fit.radial.periodic(flag)).^2./var(S.radial.hat(flag));

	obj.stat = stat;
	obj.S    = S;
end

