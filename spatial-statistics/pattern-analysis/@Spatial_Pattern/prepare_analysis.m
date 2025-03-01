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
%
%% compute quantiles of the patterns spectrum, these are used to determine
%% the target sampling resolution
function prepare_analysis(obj)
	timer = tic();

	n = obj.n;
	L = obj.L;
		
	b     = obj.b;

	if (isempty(obj.msk.b))
		obj.msk.b = true(size(b));
	else
		if (0 == max(obj.msk.b))
			error('Spatial_Pattern:EmptyArea','Masked area is empty');
		end
	end

	switch (obj.opt.datatype)
	case {'single'}
		b = single(b);
	case {'double'}
		b = double(b);
	otherwise
		error('Spatial_Pattern:DataType','Dataype must be single or double');
	end

	% convert mask to logical
	obj.msk.b = (obj.msk.b > 0);

	% resample, to make dx identical to dy
	% (depending on the grid projection, this might not be the case)
	dx = L./n;
	n_ = round(L./min(dx));
	if (n_(1)>n(1))
		% TODO implement
		% note : this is not necessary as the current input files have dx=dy	
		error('SpatialPattern','not yet implemented');
	end
	if (n_(2)>n(2))
		% TODO implement
		% note : this is not necessary as the current input files have dx=dy	
		error('Spatial_Pattern','not yet implemented');
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

	% normalize volume of the two-dimensional density to 1
	obj.S.hat = obj.S.hat/(sum(obj.S.hat,'all')*df(1)*df(2));

	% radial density (radial periodogram)
	[Sr, obj.f.r,Cr] = periodogram_radial(obj.S.hat,L);
	obj.S.radial.hat = Sr.normalized;

	% cumulative distribution
	%Cr = periodogram_cumulative(Sr.normalized,obj.f.r,'radial');
	%Cr = cumsum(cvec(obj.f.r).*cvec(obj.S.radial.hat));
	%Cr = Cr/Cr(end);

	% unique values, required for interpolation
	fdx = [true; cvec(Cr(2:end)~=Cr(1:end-1))];

	% quantiles
	obj.stat.q.fr.p05 = interp1(Cr(fdx),obj.f.r(fdx),0.05,'linear');
	obj.stat.q.fr.p50 = interp1(Cr(fdx),obj.f.r(fdx),0.50,'linear');
	obj.stat.q.fr.p84 = interp1(Cr(fdx),obj.f.r(fdx),0.84,'linear');
	obj.stat.q.fr.p95 = interp1(Cr(fdx),obj.f.r(fdx),0.95,'linear');
	obj.stat.q.fr.max = obj.f.r(end);

	obj.C.r = Cr;
	obj.b_square = b_square;
	obj.stat = setfield_deep(obj.stat,'runtime.prepare_analysis',toc(timer));
end

