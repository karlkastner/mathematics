% Sat Nov  3 12:47:26 MSK 2012
% Karl Kästner, Berlin

function [R Vr] = fem_get_2d_radial(P, T, v, n_bin, r_max)
	
	% coordinates of histrogram bins
	R = linspace(0,r_max,n_bin);
	% values of history bins
	Vr = zeros(n_bin,1);
	% weigth to norm history bin
	W  = zeros(n_bin,1);

	% get the point distancs from the origin
	Pr = sqrt(sum(P.^2,2));

	% for each triangle
	for tdx=1:size(T,1)
		% get the area
		A = [	1 P(T(tdx,1),:);
			1 P(T(tdx,2),:);
			1 P(T(tdx,3),:) ];
		area = 0.5*abs(det(A));

		% for each point of the triangle
		% add the individual components to the histogram bins
		for pdx=1:size(T,2)
			% radius
			r = Pr(T(tdx,pdx));
			% bin from radius
			b = floor(r/r_max*(n_bin-1))+1;
			% value
			Vr(b,1) = v(T(tdx,pdx),1).^2*area;
			% weight by area
			W(b,1) = W(b,1) + area;
		end % for pdx
	end % for tdx

	% unweight histogram bins by area
	Vr = Vr ./ W;
	Vr = (1./norm(Vr))*Vr;
end % fem_plot_2d_radial

