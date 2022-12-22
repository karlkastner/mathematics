% Tue  8 Nov 17:13:41 CET 2022

% determine if pattern is isotropic (spotted, labyrinthic or gapped) or
% anisotropic (banded)

function [isisotropic, resn, S] = isisotropic(Shat,L)
	[fx,fy,fr] = fourier_axis_2d(L,size(Shat));
	fxx = repmat(cvec(fx),1,length(fy));
	fyy = repmat(rvec(fy),length(fx),1);

	%for idx=1:2	
	slope = least_squares_perpendicular_offset(fxx(:),fyy(:),Shat(:));
	a = atand(slope)
	Shat = ifftshift(imrotate(fftshift(Shat),-a,'crop'));
	%end

	% reconstruct 2D density from isotropic radial density
	[Sri,fri] = periodogram_radial(Shat,[L,L]);
	S.iso     = interp1(fri,Sri.mu,fr,'linear');

	% reconstruct density from Sx and Sy
	% TODO rotate if necessary
	Sx = mean(Shat);
	Sy = mean(Shat')';
	S.aniso = Sy*Sx;
	

	% normalize
	Shat  = Shat.*L(1).*L(2)./sum(sum(Shat,1),2);
	S.iso = S.iso.*L(1).*L(2)./sum(sum(S.iso,1),2);
	S.aniso = S.aniso.*L(1).*L(2)./sum(sum(S.aniso,1),2);
	
	% residuals
	resn.iso   = rms(flat(S.iso - Shat));
	resn.aniso = rms(flat(S.aniso - Shat));
	
	isisotropic = resn.iso < resn.aniso;
end

