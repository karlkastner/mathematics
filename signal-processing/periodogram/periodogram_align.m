% Mon  5 Dec 14:28:58 CET 2022
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
function [Shat, angle_deg,out] = periodogram_align(Shat,L,mode,nf)
	if (nargin()<3||isempty(mode))
		%mode = 'max';
		mode = 'angular';
	end

	n = size(Shat);
	[fx,fy] = fourier_axis_2d(L,n);
	fxx = repmat(cvec(fx),1,n(2));
	fyy = repmat(rvec(fy),n(1),1);

	n = size(Shat);
	if (n(1)~=n(2) || abs(L(2)-L(1)) > sqrt(eps)*mean(L))
		error('periodogram must be square as it will distort by rotation')
	end

	switch (mode)
	case {'ls'}
		% least squares
		A = [ones(numel(Shat),1),fxx(:)];
		%w = Shat(:);
		W = diag(sparse(Shat(:)));
		c = (A'*W*A) \ (A'*W*fyy(:));
		slope = c(2);
		if (abs(c(2))>1)
			% swap to improve numerical accuracy
			A = [ones(numel(Shat),1),fyy(:)];
			%Shat_ = rot90(Shat);
			W = diag(sparse(Shat(:)));
			c = (A'*W*A) \ (A'*W*fxx(:));
			slope = 1./c(2);
		end
       		angle_deg = atand(slope);
	case {'po','perpendicular-offset'}
		% this somethins yields incorrectly 0
	        slope = least_squares_perpendicular_offset(fxx(:),fyy(:),Shat(:));
        	angle_deg = atand(slope);
	case {'max'}
		if (nargin()>3 && ~isempty(nf))
			Shat_ = ifftshift(trifilt2(fftshift(Shat),nf));
		else
			Shat_ = Shat;
		end
		[~,mdx] = max(Shat_,[],'all');
		
		if (fxx(mdx)==0)
			slope = inf;
		else
		slope = fyy(mdx)./fxx(mdx);
		end
        	angle_deg = atand(slope);
	case {'angular'}
		[Sa,angle] = periodogram_angular(Shat,L,nf);
		[~,mdx] = max(Sa);
		angle_rad    = angle(mdx);
		angle_deg    = rad2deg(angle_rad);
		out.Sa = Sa;
		out.angle = angle;
	otherwise
		error('unknown mode');
	end

	%n = size(Shat);
	Shat = fft_rotate(Shat,-angle_deg);
end

