% Thu Aug  9 16:22:44 MSK 2012
% Karl Kästner, Berlin

function V = vander_3d(X,p)
	for idx=1:size(X,1)
		n=1;
		Xp(1,1:3) = sym(1);
		for jdx=1:p
			Xp(jdx+1,1) = Xp(jdx,1)*X(idx,1);
			Xp(jdx+1,2) = Xp(jdx,2)*X(idx,2);
			Xp(jdx+1,3) = Xp(jdx,3)*X(idx,3);
		end

		for pdx=0:p
			for xp=pdx:-1:0
				for yp=pdx-xp:-1:0
					zp=pdx-xp-yp;
					%V(1,n) = X(idx,1)^xp * X(idx,2)^yp * X(idx,3)^zp;
					V(idx,n) = Xp(xp+1,1) * Xp(yp+1,2) * Xp(zp+1,3);
					n = n+1;
				end % yp
			end % xp
		end % pdx
	end % idx
end % vander_3d

