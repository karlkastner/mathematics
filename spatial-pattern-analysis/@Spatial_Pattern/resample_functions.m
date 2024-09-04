% 2023-01-11 14:15:11.060721307 +0100
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
%% resample empirical densities to a comman grid
%
function [Si,Ri,obj] = resample_functions(obj,xi,fi)
	% n.b. for higher than linear, posisitivty of density can be violated
	% TODO make class member
	imethod = 'linear';
	% resample 1d-densities
	% for computing the joint normalized density of.pdf patterns, the patterns
	% are resampled to a common grid interval

	if (obj.stat.isisotropic)
		lc = 1./obj.stat.fc.radial.(obj.opt.scalefield);
	else
		lc = 1./obj.stat.fc.x.(obj.opt.scalefield);
	end
	
	for field={'hat','hp','bar'}
		if (isfinite(lc))
			% isotropic
			iSr      = cumint_trapezoidal(obj.f.r,obj.S.radial.(field{1}));
			% Note that when the cdf is scaled, iSr is implicitely scaled and there is no explicit division by lc!!!
			% extrapolation is at the right hand with 1 as lim F->inf = 1
			Si.radial.cdf.(field{1}) = interp1(obj.f.r*lc,iSr,fi.radial,imethod,1);
			Si.radial.pdf.(field{1}) = derivative1(fi.radial,double(Si.radial.cdf.(field{1})));
			iSa = cumint_trapezoidal(obj.f.angle,obj.S.rot.angular.(field{1}));
			Si.angular.cdf.(field{1}) = interp1(obj.f.angle,iSa,fi.angular,'linear');
			Si.angular.pdf.(field{1}) = derivative1(fi.angular,double(Si.angular.cdf.(field{1})));
			fdai = (fi.angular>=-pi/2 & fi.angular<pi/2);
			Si.angular_p.pdf.(field{1}) = 2*cvec(Si.angular.pdf.(field{1})).*cvec(fdai);
			Ri.radial.(field{1})  = interp1(obj.r/lc,obj.R.radial.(field{1}),xi,imethod,0); 

			% anisotropic
			fdx      = (obj.f.x>=0);
			fdxi     = (fi.x>=0);
			iSx      = cumint_trapezoidal(obj.f.x(fdx), obj.S.rot.x.(field{1})(fdx));
			Si.x.cdf.(field{1}) = interp1(obj.f.x(fdx)*lc,iSx,fi.x,imethod,1);
			Si.x.pdf.(field{1}) = derivative1(fi.x,double(Si.x.cdf.(field{1})));
			Si.xp.pdf.(field{1}) = 2*Si.x.pdf.(field{1}).*(fdxi);
			n        = length(obj.R.rot.x.(field{1}));
			x_       = fourier_axis(1,n);
			Ri.x.(field{1}) = interp1(fftshift(x_)/lc,fftshift(obj.R.rot.x.(field{1})),xi,imethod,0);
			% full axis for y
			%fdx      = (obj.f.y>=0);
			%fdx      = true(size(obj.f.y));
			[void,fdx] = sort(obj.f.y);
			iSy      = cumint_trapezoidal(obj.f.y(fdx), obj.S.rot.y.(field{1})(fdx));
			Si.y.cdf.(field{1}) = interp1(obj.f.y(fdx)*lc,iSy,fi.y,imethod,1);
			Si.y.pdf.(field{1}) = derivative1(fi.y,double(Si.y.cdf.(field{1})));
		else
			Si.x.pdf.(field{1})       = NaN(1,length(fi.x));
			Si.xp.pdf.(field{1})       = NaN(1,length(fi.x));
			Si.x.cdf.(field{1})       = NaN(1,length(fi.x));
			Si.y.pdf.(field{1})       = NaN(1,length(fi.y));
			Si.y.cdf.(field{1})       = NaN(1,length(fi.y));
			Si.radial.cdf.(field{1})  = NaN(1,length(fi.radial));
			Si.radial.pdf.(field{1})  = NaN(1,length(fi.radial));
			Si.angular.cdf.(field{1}) = NaN(1,length(fi.angular));
			Si.angular.pdf.(field{1}) = NaN(1,length(fi.angular));
			Si.angular_p.pdf.(field{1}) = NaN(1,length(fi.angular));
			Ri.x.(field{1})           = NaN(1,length(xi));
			Ri.radial.(field{1})      = NaN(1,length(xi));
		end % else of if isfinite lc
		[obj.stat.Sc.angular_resampled.pdf.(field{1}),id] = max(Si.angular.pdf.(field{1}));
	        obj.stat.fc.angular_resampled.pdf.(field{1})      = fi.angular(id);
	end % for field
end % function

