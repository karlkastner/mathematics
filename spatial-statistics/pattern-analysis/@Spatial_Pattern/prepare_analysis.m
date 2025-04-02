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
	if (isempty(L))
		warning('Physical dimensions not given, using pixels');
		L = n;
	end

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
	obj.b_.square = b;
	obj.msk.b_square = obj.msk.b;

	% resize, to make domain square
	nmax = max(n);
	L    = L.*nmax./n;
 	if (n(1) < nmax)
		n(1) = nmax;
		obj.b_.square(nmax,1) = 0;
		obj.msk.b_square(nmax,1) = false;
	end
	if (n(1) > n(2))
		nmax = n(1);
		obj.b_.square(1,nmax) = 0;
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
	stat.mean_b  = wmean(double(obj.msk.b_square(:)),obj.b_.square(:));

	% subtract mean
	% note that his is superfluluous, since S(0) is later set to zero
	obj.b_.square  = (obj.b_.square - stat.mean_b);

	% rms = std, because mean has been subtracted
	stat.sd_b = wrms(obj.msk.b_square(:),obj.b_.square(:));

	% normalize
	obj.b_.square = obj.b_.square/stat.sd_b;

	% Fourier transform
	fb = fft2(obj.msk.b_square.*obj.b_.square);
	if (~isempty(obj.source))
		fe = fft2(obj.source);
	end

	% 2D periodogram, not normalized
	obj.S.hat    = conj(fb).*fb;
	if (~isempty(obj.source))
		% source spectrum
		obj.S.e.hat  = conj(fe).*fe;
		% cross-spectrum between pattern and source
		obj.S.be.hat = conj(fb).*fe; 
	end

	%if (~isempty(obj.source.S))
	%	obj.S.hat = obj.S.hat./obj.source.S;
	%	obj.S.hat(1) = 0;
	%	fdx = obj.source.S == 0;
	%	obj.S.hat(fdx) = 0;
	%end

	% normalize volume of the two-dimensional density to 1
	iSb = (sum(obj.S.hat,'all')*df(1)*df(2));
	obj.S.hat = obj.S.hat/iSb;
	if (~isempty(obj.source))
		iSe = (sum(obj.S.e.hat,'all')*df(1)*df(2));
		obj.S.e.hat = obj.S.e.hat/iSe;
		obj.S.be.hat = obj.S.be.hat/sqrt(iSb*iSe);
	end

	% radial density (radial periodogram)
	[Sbr, obj.f.r,Cr] = periodogram_radial(obj.S.hat,L);
	obj.S.radial.hat = Sbr.normalized;
	if (~isempty(obj.source))
		Ser = periodogram_radial(obj.S.e.hat,L);
		obj.S.e.radial  = Ser.normalized;
		Sber = periodogram_radial(obj.S.be.hat,L);
		% TODO normalization of the cross spectrum
		obj.S.be.radial = Sber.mu;
		% transfer function
		obj.T.radial = obj.S.be.radial./obj.S.e.radial;
		obj.S.coherence.radial = abs(Sber.mu).^2./(Sbr.mu.*Ser.mu)
	end

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
	obj.stat = setfield_deep(obj.stat,'runtime.prepare_analysis',toc(timer));
end

