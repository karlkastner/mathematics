% Tue  8 Nov 17:13:41 CET 2022
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
% TODO apply iteratively
% determine if pattern is isotropic (spotted, labyrinthic or gapped) or
% anisotropic (banded)
function [isisotropic,stat,S,f] = separate_isotropic_from_anisotropic_density(Shat,fmsk,L,varargin)
	n = size(Shat);
	kmax = 1;

	% output
	S = struct();
	f = struct();

	% apply mask
	S.hat = Shat.*fmsk;
	% normalize
	S.hat   = S.hat.*L(1).*L(2)./sum(sum(S.hat,1),2);

	% rotate coordinate system so that maximum occurs on x-axis
	[S.hat,stat.angle_deg] = periodogram_align(Shat,L,varargin{:});

	% mask extra-diagonal pattern
	%x = zeros(100,100);
	%[fx,fy,fr] = fourier_axis_2d([1,1],[100,100]);
	%S.msk = (abs(f.x)>abs(f.y'));
	%S.hat_iso_msk = S.hat.*(~S.msk);

	% initial guess
	S.hat_iso   = S.hat;
	S.hat_aniso = S.hat;
	%p = 0;

	% fixed point iteration
	for kdx=1:kmax

	% estimate radial density
	[S.r,f.r,~,A] = periodogram_radial(S.hat_iso,L);

	% reconstruct 2D density from isotropic radial density
	S.iso = reshape(A*S.r.normalized,n);

	% anisotropic densities
	S.x = mean(S.hat_aniso,2);
	S.y = mean(S.hat_aniso,1)';
%	S.x = S.x/sum(S.x(f.x>=0)*(f.x(2)-f.x(1)));
%	S.y = S.y/sum(S.y(f.y>=0)*(f.y(2)-f.y(1)));
	% reconstruct 2D density from Sx and Sy
	S.aniso = S.x*S.y';

	% normalize
	S.iso   = S.iso.*L(1).*L(2)./sum(sum(S.iso,1),2);
	S.aniso = S.aniso.*L(1).*L(2)./sum(sum(S.aniso,1),2);

	% normalized residuals
	stat.resn.iso   = rms(S.iso  - S.hat,'all');
	stat.resn.aniso = rms(S.aniso - S.hat,'all');

	% this is not a p value and in range (-1,1)
	% 100% aniso p = (1+(0-1))/2 = 0
	% 100% iso   p = (1+(1-0))/2 = 1
	%stat.piso   = (1 + (stat.resn_aniso - stat.resn_iso)./rms(S.hat,'all'))/2.0;
	
	% Least squares:
	% S_hat = p*S_iso + (1-p)*S_aniso
	% (S_iso - S_aniso) p = (S_hat - S_aniso)
	stat.p_iso = flat(S.iso - S.aniso) \ flat(S.hat - S.aniso);
	isisotropic = stat.p_iso>0.5;


%	% factor for anisotropic part
%	% alternatively, we can raise S to the power
%	if (0)
%	fa = 2.5; % for 2, aniso are missclassified as iso, for 3, iso are misscl. as ansio
%	piso = (sum(S.iso(:)) - fa*sum(S.hat_aniso(:)))./sum(Shat(:));
%	else
%		p = 1.5;
%		e_iso   = 2*sum(flat(S.hat_iso_msk.^p));
%		e_aniso =   sum(flat(max(S.hat_aniso,0).^p));
%		piso = (e_iso - e_aniso)./(e_iso + e_aniso);
%	end

	%S.hat_iso   = Shat - (1-p)*S.aniso*rms(S.hat,'all')./rms(S.aniso,'all');
	%S.hat_aniso = Shat - p*S.iso*rms(S.hat,'all')./rms(S.iso,'all');

	% anisotropic and isotropic part of ther periodogram,
	% note : this can yield negative values
	S.hat_iso   = Shat - (1-stat.p_iso)*S.aniso;
	S.hat_aniso = Shat - stat.p_iso*S.iso;

	end % for kdx
end % function

