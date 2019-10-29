% Thu Aug  9 16:25:12 MSK 2012
% Karl Kästner, Berlin

function [C_dx C_dy C_dz] = derivative_3d(C,p)
	for idx=1:size(C,1)
		nx=1;
		ny=1;
		nz=1;
		n=1;
		for pdx=0:p
			for xp=pdx:-1:0
				for yp=pdx-xp:-1:0
					zp=pdx-xp-yp;
					if (xp > 0)
						C_dx(idx,nx) = xp*C(idx,n);
						nx = nx+1;
					end
					if (yp > 0)
						C_dy(idx,ny) = yp*C(idx,n);
						ny = ny+1;
					end
					if (zp > 0)
						C_dz(idx,nz) = zp*C(idx,n);
						nz = nz+1;
					end
					n = n+1;
				end % yp
			end % xp
		end % pdx
	end % for idx
end % derivative_3d

