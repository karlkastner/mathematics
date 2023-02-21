% 2021-06-21 22:52:10.758702394 +0200
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution
%
% Sxy : nxn
% Sra : n/2*(pi*n/2)
function [Sra,fr,angle,out,A] = fft2_cartesian2radial(Sxy,L)
	nxy = size(Sxy);
	% TODO allow for asymmetric matrices
	% requires appropriate scaling of the coodinates
	if (nxy(2) ~= nxy(1) || abs(L(1)-L(2))>sqrt(L(1)*L(2)*eps))
		error('matrix must be square');
	end

	% size of Sra
	% discretizing the maximum circumference ensures that sectors
	% at the maximum radius are only 1 pixel wide
	nra = [floor(nxy(1)/2)-1, round(pi*nxy(1))];
	% radii, as a fraction of L/2
	rid = (1:nra(1))';
	% radial frequencies
	dfr = 1./L(1);
	fr  = (rid-1)*dfr;
	% angle, as a fraction of 2 pi
	aid = 1:nra(2);
	% angle
	da    = 2*pi/nra(2);
	angle = da*(aid-1)-pi;
	% sine of angle
	s = sin(angle);
	% cosine of angle
	c = cos(angle);
	% x-index
	x = rid*c;
	y = rid*s;
	% integer parts
	% note, this really has to be floor not fix
	i = floor(x);
	j = floor(y);
	% fractional parts for bilinear interpolation
	p = 1-(x-i);
	q = 1-(y-j);

	% bilinear interpolation
	Sra = (  Sxy(sub2ind(nxy,fourier_freq2ind(i,nxy(1)),fourier_freq2ind(j,nxy(2)))).*p.*q ...
                   + Sxy(sub2ind(nxy,fourier_freq2ind(i+1,nxy(1)),fourier_freq2ind(j,nxy(2)))).*(1-p).*q ...
		   + Sxy(sub2ind(nxy,fourier_freq2ind(i,nxy(1)),fourier_freq2ind(j+1,nxy(2)))).*p.*(1-q) ...
		   + Sxy(sub2ind(nxy,fourier_freq2ind(i+1,nxy(1)),fourier_freq2ind(j+1,nxy(2)))).*(1-p).*(1-q) ...
	        );

	out.i = i;
	out.j = j;
	out.p = p;
	out.q = q;
	out.nra = nra;
	out.dfr = dfr;
	out.da  = da;

	% normalize
	%w   = fr./(sum(fr)*dfr);
	

	% TODO matrix output
	if (nargout()>1)
	end
end

