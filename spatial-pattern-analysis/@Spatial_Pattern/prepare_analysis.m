% Wed 18 May 13:50:47 CEST 2022
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
function prepare_analysis(obj)
	n = obj.n;
	L = obj.L;
		
	b     = obj.b;
	%msk.b = obj.msk.b;
	if (isempty(obj.msk.b))
		obj.msk.b = true(size(b));
	else
		if (0 == max(obj.msk.b))
			error('Spatial_Pattern','Masked area is empty');
		end
	end

	% convert image to double
	b    = double(b);

	% convert mask to logical
	obj.msk.b = (obj.msk.b > 0);

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

	% make square
	b_square = b;
	obj.msk.b_square = obj.msk.b;

	% resize, to make domain square
	nmax = max(n);
	L    = L.*nmax./n;
 	if (n(1) < nmax)
		n(1) = nmax;
		b_square(nmax,1) = 0;
		obj.msk.b_square(nmax,1) = false;
	end
	if (n(1) > n(2))
		nmax = n(1);
		b_square(1,nmax) = 0;
		obj.msk.b_square(1,nmax) = false;
	end
	df   = 1./L;
	n    = [nmax,nmax];
	obj.stat.L_square = L;
	obj.stat.n_square = n;
	df = 1./obj.stat.L_square;

	% grid in real space
	% note : x an y do not match original figure when the figure was resampled or or resized
	obj.x   = linspace(-L(1)/2,L(1)/2,n(1))';
	obj.y   = linspace(-L(2)/2,L(2)/2,n(2))';

	% grid in frequency space
	[obj.f.x,obj.f.y,obj.f.rr,obj.f.tt] = fourier_axis_2d(L,n);

	% the weighted mean allows for feathering the transition
	stat.mean_b  = wmean(double(obj.msk.b_square(:)),b_square(:));

	% subtract mean
	b_square  = (b_square - stat.mean_b);

	% rms = std, because mean has been subtracted
	stat.sd_b = wrms(obj.msk.b_square(:),b_square(:));

	% normalize
	b_square = b_square/stat.sd_b;

	% 2D periodogram, not yet normalized
	obj.S.hat  = abs(fft2(obj.msk.b_square.*b_square)).^2;
	% normalize
	obj.S.hat = obj.S.hat/(0.5*sum(obj.S.hat,'all')*df(1)*df(2));

	% radial density (radial periodogram)
	[Sr, obj.f.r] = periodogram_radial(obj.S.hat,L);
	obj.S.radial.hat = Sr.normalized;

	% cumulative distribution
	Cr = cumsum(cvec(obj.f.r).*cvec(obj.S.radial.hat));
	Cr = Cr/Cr(end);

	% make values unique (quick hack)
	Cr = cvec(Cr) + (0:length(Cr)-1)'*1e-12; 

	% quantiles
	obj.stat.q.fr.p05 = interp1(Cr,obj.f.r,0.05,'linear');
	obj.stat.q.fr.p50 = interp1(Cr,obj.f.r,0.50,'linear');
	obj.stat.q.fr.p84 = interp1(Cr,obj.f.r,0.84,'linear');
	obj.stat.q.fr.p95 = interp1(Cr,obj.f.r,0.95,'linear');
	obj.stat.q.fr.max = obj.f.r(end);

	obj.C.r = Cr;
	obj.b_square = b_square;
end

