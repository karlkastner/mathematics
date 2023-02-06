% 2023-01-11 14:15:11.060721307 +0100
% TODO:
%	compute cdf (with fourier?)
%	Fourier interpolation of cdf
%	restore monotonicity by sorting
%	compute pdf as Fourier derivative -> can this still result in negative values?
% non-fourier:
%	compute integral as trapezoidal rule (1/2 U + 1/2L) and derivative with central differences
%	as mid-point and lower/upper induce a phase shift!!
function [Si,Ri,obj] = resample_functions(obj,xi,fi)
	% n.b. for higher than linear, posisitivty of density can be violated
	imethod = 'linear';
	% resample 1d-densities
	% for computing the joint normalized density of.pdf patterns, the patterns
	% are resampled to a common grid interval
	Si.x.pdf = NaN(1,length(fi.x));
	Si.x.cdf = NaN(1,length(fi.x));
	Si.y.pdf = NaN(1,length(fi.y));
	Si.y.cdf = NaN(1,length(fi.y));
	Si.radial.cdf = NaN(1,length(fi.x));
	Si.radial.pdf = NaN(1,length(fi.x));
	Si.angular.cdf = NaN(1,length(fi.angular));
	Si.angular.pdf = NaN(1,length(fi.angular));
	Ri.x          = NaN(1,length(xi));
	Ri.radial     = NaN(1,length(xi));

	if (obj.stat.isisotropic)
		lc = 1./obj.stat.fc.radial.clip;
	else
		lc = 1./obj.stat.fc.x.clip;
	end
	
% TODO repeat for hat and bar
	if (isfinite(lc))
		% isotropic
		iSr      = cumint_trapezoidal(obj.f.r,obj.S.radial.clip);
		% Note that when the cdf is scaled, iSr is implicitely scaled and there is no explicit division by lc!!!
		% extrapolation is at the right hand with 1 as lim F->inf = 1
		Si.radial.cdf = interp1(obj.f.r*lc,iSr,fi.x,imethod,1);
		Si.radial.pdf = derivative1(fi.x,Si.radial.cdf);
		%iSa = [0;cumsum([obj.S.angular.clip; obj.S.angular.clip(1)])];
		%Si.angular = interp1(f,S,fi.angular,'linear');
		iSa = cumint_trapezoidal(obj.f.angle,obj.S.angular.clip);
		Si.angular.cdf = interp1(obj.f.angle,iSa,fi.angular,'linear');
		Si.angular.pdf = derivative1(fi.angular,Si.angular.cdf);
		%Si.angular.pdf = diff(Si.angular.cdf)./diff(cvec(fi.angle));
		%Si.angular.pdf = 2*Si.angular.pdf./sum(mid(Si.angular).*diff(fi.angle));
	
		Ri.radial  = interp1(obj.r/lc,obj.R.radial.clip,xi,imethod,0); 

		%  obj.S.angular.clip;
		% anisotropic
		%fx       = obj.f.x;
		fdx      = (obj.f.x>=0);
		iSx      = cumint_trapezoidal(obj.f.x(fdx), obj.S.rot.x.clip(fdx));
		%[0;cumsum(obj.S.rot.x.clip(fdx))*(fx(2)-fx(1))];
		%Si.x.cdf = interp1(inner2outer(fx(fdx))*lc,iSx/lc,inner2outer(fi.x),imethod,0);
		Si.x.cdf = interp1(obj.f.x(fdx)*lc,iSx,fi.x,imethod,1);
		Si.x.pdf = derivative1(fi.x,Si.x.cdf);
		%./(fi.x(2)-fi.x(1));
		n        = length(obj.R.rot.x.clip);
		x_       = fourier_axis(1,n);
		Ri.x     = interp1(fftshift(x_)/lc,fftshift(obj.R.rot.x.clip),xi,imethod,0);
		fdx      = (obj.f.y>=0);
		iSy      = cumint_trapezoidal(obj.f.y(fdx), obj.S.rot.y.clip(fdx));
		%yi      = linspace(0,3.5,length(obj.S.rot.y.clip));
		Si.y.cdf = interp1(obj.f.y(fdx)*lc,iSy,fi.y,imethod,1);
		Si.y.pdf = derivative1(fi.y,Si.y.cdf);

	end % if isfinite lc
end % function

