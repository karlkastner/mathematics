% Thu May 24 19:53:08 MSK 2012
% Karl Kästner, Berlin

function [Mx My err_max err] = fdm_mark_unstructured_2d(P, v)
	np = size(P,1);
	Mx = zeros(np,2);
	My = zeros(np,2);
	nx = 0;
	ny = 0;

	% threshhold for refinement
	p = 0.5;

	% first order partial derivative matrices
	D1x = fdm_assemble_d1_2d(P, N, [1 0]);
	D1y = fdm_assemble_d1_2d(P, N, [0 1]);
	% second order partial derivatives matrices
	D2x = fdm_assemble_d2_2d(P, N, [1 0]);
	D2y = fdm_assemble_d2_2d(P, N, [0 1]);
	% second order partial derivatives
	v2x = D2x*v;
	v2y = D2y*v;
	% third order partial derivatives
	v3x = D1x*v2x;
	v3y = D1x*v2y;
	% fourth order partial derivatives
	v4x = D2x*v2x;
	v4y = D2y*v2y;

	% step widths
	[Hxl Hxr Hyl Hyr] = fdm_h_unstructured_2d(P);

	% approximated error
	errx = 1/3*abs((Hxr - Hxl).*v3x) + 1/12*abs((Hxr.^2 - Hxl.*Hxr + Hxr.^2).*v4x));
	erry = 1/3*abs((Hyr - Hyl).*v3y) + 1/12*abs((Hyr.^2 - Hyl.*Hyr + Hyr.^2).*v4y));

	% maximum error (for splitting criteria)
	err_max = max(max(errx), max(erry));

	% total error independend of axis (for plot)
	if (nargout() > 3)
		err = errx+erry;
	end

% TODO L1 not inf norm, otherwise no convergence
	% select points for refinement
	for idx=1:np
		% check gap to right neighbour
		% no need to check gap to left neighbour, will be done by left neighbour
		if (0 == N(idx,6) && N(idx,2) > 0)
			% right neighbour is unique and exists
			if ( 0.5*(errx(idx,1) + errx(N(idx,2),1)) > p*err_max)
				% left and right error of line segment is high -> split
				nx = nx+1;
				Mx(nx,:) = [idx, N(idx,2)];
			end
		else

		% top neighbour
		if (0 == N(idx,8) && N(idx,4) > 0)
			% top neighbour is unique and exists
			if ( 0.5*(errx(idx,1) + errx(N(idx,4),1)) > p*err_max)
				ny = ny+1;
				My(ny,:) = [idx, N(idx,4)];
			end
		else
	end % for idx (all points)

	% trim to actual size
	Mx = Mx(1:nx,:);
	My = My(1:ny,:);
end % fdm_mark_unstructured_2d()

