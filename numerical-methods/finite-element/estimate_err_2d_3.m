% Thu Jun 14 16:59:36 MSK 2012
% Karl Kästner, Berlin

function [err nH] = estimate_err_2d_3(selection, N, dV, area, h_side, s_angle, C, err, nH)
	for tdx=selection
		% calculate seminorm of second derivative
		for jdx=1:3
			a = tdx;
			b = N(tdx,jdx);
			% check that neighbour exists (not a domain boundary)
			if (b > 0)
				dx      = C(a,1) - C(b,1);
				dy      = C(a,2) - C(b,2);
				dr_sqr  = dx^2 + dy^2;
				dv_x    = dV(a,1) - dV(b,1);
				dv_y    = dV(a,2) - dV(b,2);
				dv_xx   = dv_x*dx/dr_sqr;
				dv_xy   = dv_x*dy/dr_sqr;
				dv_yx   = dv_y*dx/dr_sqr;
				dv_yy   = dv_y*dy/dr_sqr;
				% Hessian matrix
				H = [            dv_xx  0.5*(dv_xy + dv_yx);
                                      0.5*(dv_xy + dv_yx)            dv_yy ];
				% derivative norm
				nH(tdx) = max(nH(tdx), norm(H,'inf'));
			end % if b>0
		end % jdx
		% error contribution
		err(tdx) = 1/6 * nH(tdx) * (prod(h_side(tdx,:))/prod(s_angle(tdx,:)))^(2/3);
		%err(tdx) = 1/6 * nH(tdx) * (max(h_side(tdx,:))/min(s_angle(tdx,:)))^2;
	end % for tdx
end % estimate_err_2d_3

