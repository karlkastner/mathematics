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
%% input:
%% 	Sxy : nxn
%% ouput:
%% 	Sra : n/2*(pi*n/2)
%%	angle
%% function [Sa,angle,A] = periodogram_angular(Sxy,L,nf)
function [Sa,angle,A] = periodogram_angular(Sxy,L,nf)

	if (nargin()<3)
		nf = [];
	end

	[Sra,fr,angle,out] = fft2_cartesian2radial(Sxy,L);
	

	% density estimate, smooth along angles
	if (~isempty(nf))
		Sra = gaussfilt1(Sra',nf)';
	end

	% the low frequency bins are part of several sectors and thus are weighted down in proportion
	Sa = Sra'*fr;

	% matrix output
	if (nargout()>2)
	nxy = size(Sxy);
	%nra = size(Sra);
	i = out.i;
	j = out.j;
	p = out.p;
	q = out.q;
	nra = out.nra;

	buf = [ i(:),j(:)  ,flat(p.*q.*fr);
                i(:)+1,j(:)  ,flat((1-p).*q.*fr);
		i(:)  ,j(:)+1,flat(p.*(1-q).*fr);
                i(:)+1,j(:)+1,flat((1-p).*(1-q).*fr)
	      ];
	ind_Sa   = repmat(flat(repmat(1:nra(2),nra(1),1)),4,1);

	% translate indices in the same manner like ifftshift
	buf(:,1) = fourier_freq2ind(buf(:,1),nxy(1));
	buf(:,2) = fourier_freq2ind(buf(:,2),nxy(2));

	% indices in flattened matrix
	ind_Sxy  = sub2ind(nxy,buf(:,1),buf(:,2));	
	% averaging matrix
	A        = sparse(ind_Sa,ind_Sxy,buf(:,3),nra(2),nxy(1)*nxy(2));
if (0)
	% normalize, note this is up to round-off already done
	% TODO normalize columns
	s = sum(A,2);
	A = A .*(1./(pi*s));
end
	% average spectrum
	if (0)
		Sa = A*flat(Sxy);
	end
	end
	% normalize spectrum over the full-circle
	Sa = Sa./(sum(Sa)*out.da);
	
end

