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

	% without mean component
	[par0(1),par0(2)] = spectral_density_brownian_phase_mode2par(fc.x.bar,Sc.x.bar);
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,'brownian-phase',par0,obj.opt.objective,nf);
	Sfit = Sfit/sum(Sfit)*sum(S.rot.x.hat);
	S.rot.x.brownian_phase = Sfit;
	stat.fit.x.brownian_phase.par  = par;
	stat.fit.x.brownian_phase.stat = fitstat;

	% including mean component
	par0 = [par, 5./fc.rr.bar];
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,'brownian-phase-mean',par0,obj.opt.objective,nf);
	Sfit = Sfit/sum(Sfit)*sum(S.rot.x.hat);
	S.rot.x.brownian_phase_mean = Sfit;
	stat.fit.x.brownian_phase_mean.par  = par;
	stat.fit.x.brownian_phase_mean.stat = fitstat;

	% without mean component
	par0 = [fc.rr.bar, spectral_density_bandpass_continuous_max2par(fc.x.bar,Sc.x.bar,10)];
	%par = [fc.rr.bar,20];
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,'bandpass-continuous',par0,obj.opt.objective,nf);
	Sfit = Sfit/sum(Sfit)*sum(S.rot.x.hat);
	S.rot.x.bandpass         = Sfit;
	stat.fit.x.bandpass.par  = par;
	stat.fit.x.bandpass.stat = fitstat;

	% with mean
	par0 = [par,10./fc.rr.bar];
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.x.hat,obj.w.x,'bandpass-continuous-mean',par0,obj.opt.objective,nf);
	Sfit = Sfit/sum(Sfit)*sum(S.rot.x.hat);
	S.rot.x.bandpass_mean = Sfit;
	stat.fit.x.bandpass_mean.par  = par;
	stat.fit.x.bandpass_mean.stat = fitstat;

	% fit parametric density models in the direction perpendicular to the bands
	fa = obj.f.y;
	par0 = spectral_density_brownian_phase_across_mode2par(S.rot.y.hat(1));
	[pary,Sfit,fitstat] = fit_spectral_density(obj.f.y,S.rot.y.hat,obj.w.y,'brownian-phase-across',par0,obj.opt.objective,nf);
	S.rot.y.brownian_phase_across = Sfit;
	stat.fit.y.brownian_phase_across.par  = pary;
	stat.fit.y.brownian_phase_across.stat = fitstat;
	%aS_ = aS_/(sum(aS_)*(fa(2)-fa(1)));

	% fit parametric density models to the radial density
	par0 = [];
	[par0(1),par0(2)] = spectral_density_brownian_phase_mode2par(fc.radial.bar,Sc.radial.bar);
	%Lr = 1./(obj.f.r(2)-obj.f.r(1));
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,'brownian-phase',par0,obj.opt.objective,nf);
	%cS_ = cS_/sum(cS_)*sum(cS);
	S.radial.brownian_phase = Sfit;
	stat.fit.radial.brownian_phase.par  = par;
	stat.fit.radial.brownian_phase.stat = fitstat;

	par0 = [par0,5./fc.rr.bar];
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,'brownian-phase-mean',par0,obj.opt.objective,nf);
	%cS_ = cS_/sum(cS_)*sum(cS);
	S.radial.brownian_phase_mean = Sfit;
	stat.fit.radial.brownian_phase_mean.par  = par;
	stat.fit.radial.brownian_phase_mean.stat = fitstat;

	%par = [fc.rr.bar,10];
	par0 = [fc.radial.bar, spectral_density_bandpass_continuous_max2par(fc.radial.bar,Sc.radial.bar,10)];
	%nf = 1;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,'bandpass-continuous',par0,obj.opt.objective,nf);
	S.radial.bandpass = Sfit;
	stat.fit.radial.bandpass.par  = par;
	stat.fit.radial.bandpass.stat = fitstat;

	par0 = [par,10./fc.rr.bar];
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.hat,obj.w.r,'bandpass-continuous-mean',par0,obj.opt.objective,nf);
	S.radial.bandpass_mean = Sfit;
	stat.fit.radial.bandpass_mean.par  = par;
	stat.fit.radial.bandpass_mean.stat = fitstat;

	obj.stat = stat;
	obj.S    = S;
end

