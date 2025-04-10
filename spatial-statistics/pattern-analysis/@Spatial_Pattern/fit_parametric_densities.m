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
%
% TODO the phase drift and bandpass are already "mirrored" but should be normalized to 1, not 2
%      note that this does not effect the fit as only the positive axis is chosen for the fit
% TODO add weibull and generalized gamma (3 parameter) to fit
function fit_parametric_densities(obj)
	S    = obj.S;
	stat = obj.stat;
	fc   = obj.stat.fc;
	Sc   = obj.stat.Sc;

	% TODO from obj
	nf = [];

	% fit to raw periodogram (axially averaged)
	fitfield = 'hat';
	% get initial paramters from consitently averaged density
	p0field  = 'con';

	% fit parametric density models in direction perpendicular to bands
	% this is computed anyway for all patterns, but only meaningfull for banded patterns

	% phase drift
%	try
	[par0(1),par0(2)] = phase_drift_pdf_mode2par(fc.x.(p0field),Sc.xp.(p0field));
	Sfun = @phase_drift_pdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.xp.(fitfield),obj.w.xp,Sfun,par0,obj.opt.objective,nf,[0,0]);
	S.rot.xp.phase_drift = Sfit;
	stat.fit.xp.phase_drift.par  = par;
	stat.fit.xp.phase_drift.stat = fitstat;
	[stat.fc.xp.phase_drift,stat.Sc.xp.phase_drift] = phase_drift_pdf_mode(par(1),par(2));
%	catch e
%		disp(e);
%		S.rot.x.phase_drift = NaN(size(S.rot.x.(fitfield)));
%		stat.fit.x.phase_drift.par = [NaN,NaN];
%		stat.fit.x.phase_drift.stat = struct();
%	end

	% bandpass
%	try
	[par0(1),par0(2)] = bandpass1dpdf_mode2par(fc.x.(p0field),Sc.xp.(p0field));
	%par0 = [fc.rr.(p0field), bandpass1dpdf_max2par(fc.x.(p0field),Sc.x.(p0field),10)];
	Sfun = @bandpass1dpdf;
	[par,Sfit,fitstat]       = fit_spectral_density(obj.f.x,S.rot.xp.(fitfield),obj.w.xp,Sfun,par0,obj.opt.objective,nf,[0,0.5]);
	S.rot.xp.bandpass         = Sfit;
	stat.fit.xp.bandpass.par  = par;
	stat.fit.xp.bandpass.stat = fitstat;
	stat.fc.xp.bandpass      = par(1);
	stat.Sc.xp.bandpass      = bandpass1dpdf_max(par(1),par(2));
%	catch e
%		disp(e);
%		S.rot.x.bandpass = NaN(size(S.rot.x.(fitfield)));
%		stat.fit.x.bandpass.par = [NaN,NaN];
%		stat.fit.x.bandpass.stat = struct();
%	end

	% log-normal
%	try
	[par0(1),par0(2)] = lognpdf_mode2par(fc.x.(p0field),Sc.xp.(p0field));
	Sfun = @lognpdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.xp.(fitfield),obj.w.xp,Sfun,par0,obj.opt.objective,nf,[-inf,0]);
	S.rot.xp.logn  = Sfit;
	stat.fit.xp.logn.par  = par;
	stat.fit.xp.logn.stat = fitstat;
	[stat.fc.xp.logn,stat.Sc.xp.logn] = lognpdf_mode(par(1),par(2));
%	catch e
%		disp(e);
%		S.rot.x.logn         = NaN(size(S.rot.x.(fitfield)));
%		stat.fit.x.logn.par  = [NaN,NaN];
%		stat.fit.x.logn.stat = struct();
%	end

	% gamma
	[par0(1),par0(2)] = gampdf_mode2par(fc.x.hp,Sc.xp.hp);
	par0(1) = max(par0(1),1);
	Sfun = @gampdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.xp.(fitfield),obj.w.xp,Sfun,par0,obj.opt.objective,nf);
	S.rot.xp.gamma  = Sfit;
	stat.fit.xp.gamma.par  = par;
	stat.fit.xp.gamma.stat = fitstat;
	[stat.fc.xp.gamma,stat.Sc.xp.gamma] = gampdf_mode(par(1),par(2));

	% mirrored-normal
	% note that for the log-normal and gamma densities, there is beside the
	% scale no difference between the mirrored and unimodal densities,
	% as the tail is not reaching into the negative half-axis
	[par0(1),par0(2)] = normalmirroredpdf_mode2par(fc.x.hp,0.5*Sc.xp.hp);
	%par0(1) = max(par0(1),1);
	Sfun = @(fx,a,b) 2*normalmirroredpdf(fx,a,b);
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.x,S.rot.xp.(fitfield),obj.w.x,Sfun,par0,obj.opt.objective,nf);
	S.rot.xp.normalmirrored  = Sfit;
	stat.fit.xp.normalmirrored.par  = par;
	stat.fit.xp.normalmirrored.stat = fitstat;
	[stat.fc.xp.normalmirrored,Sc_] = normalmirroredpdf_mode(par(1),par(2));
	stat.Sc.xp.normalmirrored = 2*Sc_;

	% TODO use three parametric densities allowing for skewness, like generalized gamma

	% fit parametric density models in the direction perpendicular to the bands
	fy = obj.f.y;
	fdx = (fy>=0);
	par0 = phase_drift_parallel_pdf_mode2par(S.rot.y.(fitfield)(1));
	Sfun = @phase_drift_parallel_pdf;
	[pary,Sfit,fitstat] = fit_spectral_density(obj.f.y,S.rot.y.(fitfield),obj.w.y,Sfun,par0,obj.opt.objective,nf);
	S.rot.y.phase_drift_parallel = Sfit;
	stat.fit.y.phase_drift_parallel.par  = pary;
	stat.fit.y.phase_drift_parallel.stat = fitstat;

	% flat density of white noise
	S.fit.x.white = mean(S.rot.x.(fitfield))*ones(size(S.rot.x.(fitfield)));
	dfx = obj.f.x(2)-obj.f.x(1);
	stat.fit.x.white.stat.goodness.hd = hellinger_distance(S.rot.x.(fitfield),S.fit.x.white,dfx,obj.w.x);
	stat.fit.x.white.stat.goodness.r2 = 1 - 2*stat.fit.x.white.stat.goodness.hd;

	% periodic
	% (one frequency component in single bin without sub or super-harmonics)
	% does not account for cases where the characteristic frequency spreads over two bins
	flag = obj.f.x>=0;
	[mv,mdx] = max(S.rot.x.(fitfield).*flag);
	S.fit.x.periodic = zeros(size(S.rot.x.(fitfield)));
	S.fit.x.periodic(mdx) = sum(S.rot.x.(fitfield).*flag);
	stat.fit.x.periodic.stat.goodness.hd = hellinger_distance(S.rot.x.(fitfield)(flag),S.fit.x.periodic(flag),dfx,obj.w.x(flag));
	stat.fit.x.periodic.stat.goodness.r2 = 1 - 2*stat.fit.x.periodic.stat.goodness.hd;

	%
	% fit parametric density models to the radial density
	%

	% phase drift
	par0 = [];
	fc_ = max(0.5*obj.f.r(2),fc.radial.(p0field));
	[par0(1),par0(2)] = phase_drift_pdf_mode2par(fc_,Sc.radial.(p0field));
	Sfun = @phase_drift_pdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.(fitfield),obj.w.r,Sfun,par0,obj.opt.objective,nf);
	S.fit.radial.phase_drift = Sfit;
	stat.fit.radial.phase_drift.par  = par;
	stat.fit.radial.phase_drift.stat = fitstat;
	[stat.fc.radial.phase_drift,stat.Sc.radial.phase_drift]= phase_drift_pdf_mode(par(1),par(2));

	% bandpass
	%par0 = [fc_, bandpass1dpdf_max2par(fc_,Sc.radial.(p0field),10)];
	[par0] = [fc.radial.(p0field),bandpass1dpdf_max2par(fc.radial.(p0field),Sc.radial.(p0field))];
	Sfun = @bandpass1dpdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.(fitfield),obj.w.r,Sfun,par0,obj.opt.objective,nf,[0,0.5]);
	S.fit.radial.bandpass = Sfit;
	stat.fit.radial.bandpass.par  = par;
	stat.fit.radial.bandpass.stat = fitstat;
	stat.fc.radial.bandpass = par(1);
	%,stat.Sc.radial.bandpass]= bandpass1dpdf_mode(par(1),par(2));
	[stat.Sc.radial.bandpass]= bandpass1dpdf_max(par(1),par(2));

	% log-normal
	try
	[par0(1),par0(2)] = lognpdf_mode2par(fc.radial.(p0field),Sc.radial.(p0field));
	Sfun = @lognpdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.(fitfield),obj.w.r,Sfun,par0,obj.opt.objective,nf,[-inf,0]);
	S.fit.radial.logn  = Sfit;
	stat.fit.radial.logn.par  = par;
	stat.fit.radial.logn.stat = fitstat;
	catch e
		S.radial.logn             = NaN(size(S.radial.(fitfield)));
		stat.fit.radial.logn.par  = [NaN,NaN];
		stat.fit.radial.logn.stat = struct();
	end
	[stat.fc.radial.logn,stat.Sc.radial.logn]= lognpdf_mode(par(1),par(2));

	% gamma
	[par0(1),par0(2)] = gampdf_mode2par(fc_,Sc.radial.(p0field));
	Sfun = @gampdf;
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r(2:end),S.radial.(fitfield)(2:end),obj.w.r(2:end),Sfun,par0,obj.opt.objective,nf);
	S.fit.radial.gamma  = Sfit;
	stat.fit.radial.gamma.par  = par;
	stat.fit.radial.gamma.stat = fitstat;
	[stat.fc.radial.gamma,stat.Sc.radial.gamma]= gampdf_mode(par(1),par(2));

	% mirrored normal
	% the factor 1/2 accounts for that only the positive half is fitted
	[par0(1),par0(2)] = normalmirroredpdf_mode2par(fc.radial.(p0field),0.5*Sc.radial.(p0field));
	Sfun = @(fr,a,b) 2*normalmirroredpdf(fr,a,b);
	[par,Sfit,fitstat] = fit_spectral_density(obj.f.r,S.radial.(fitfield),obj.w.r,Sfun,par0,obj.opt.objective,nf);
	S.fit.radial.normalmirroredpdf  = Sfit;
	stat.fit.radial.normalmirrored.par  = par;
	stat.fit.radial.normalmirrored.stat = fitstat;
	[stat.fc.radial.normalmirrored,Sc_]= normalmirroredpdf_mode(par(1),par(2));
	stat.Sc.radial.normalmirrored = 2*Sc_;

	% flat density of white noise
	S.fit.radial.white = mean(S.radial.(fitfield))*ones(size(S.radial.(fitfield)));
	dfr = obj.f.r(2)-obj.f.r(1);
	stat.fit.radial.white.stat.goodness.hd = hellinger_distance(S.radial.(fitfield),S.fit.radial.white,dfr,obj.w.r);
	stat.fit.radial.white.stat.goodness.r2 = 1 - 2*stat.fit.radial.white.stat.goodness.hd;

	% periodic
	flag = obj.f.r>=0;
	[mv,mdx] = max(S.radial.(fitfield).*flag);
	S.fit.radial.periodic = zeros(size(S.radial.(fitfield)));
	S.fit.radial.periodic(mdx) = sum(S.radial.(fitfield).*flag);
	stat.fit.radial.periodic.stat.goodness.hd = hellinger_distance(S.radial.(fitfield)(flag),S.fit.radial.periodic(flag),dfr,obj.w.r(flag));
	stat.fit.radial.periodic.stat.goodness.r2 = 1 - 2*stat.fit.radial.periodic.stat.goodness.hd;

	obj.stat = stat;
	obj.S    = S;
end % function fit_parametric_densities


