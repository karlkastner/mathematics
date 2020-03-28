% Thu Jun 14 22:00:47 MSK 2012
% 
% Brown's lense cusshion distrotion correction
% pin cussion distortion
function P_ = distort(P,K_r,K_t,L,mode,flag,p)
	rs = sum(P.^2,2);
	r  = sqrt(rs);
	r_max = sqrt(sum((L/2).^2));
	switch (mode)
	case {0}
		% power of r by keeping direction
		P_ = [ (P(:,1)./(eps+r)).*(r/r_max).^p ...
		       (P(:,2)./(eps+r)).*(r/r_max).^p ];
	case{1}
		% Brown's distortion
		P_ = P.*(1 + (K_r(1)*rs + K_r(2)*rs.^2)*ones(1,size(P,2))) ...
			+ [K_t(1)*(rs + 2*P(:,1).^2) + 2*K_t(2)*P(:,1).*P(:,2), ...
			   K_t(2)*(rs + 2*P(:,2).^2) + 2*K_t(1)*P(:,1).*P(:,2)];
		P_ = P_/(1 + K_r(1)*r_max^2 + K_r(2)*r_max^4);
	end
	if (flag)
		% stretch grid only in centre, keep points on boundary at old positions
		w = (max([(2*abs(P(:,1)))/L(1) (2*abs(P(:,2)))/L(2)],[],2)*ones(1,2)).^p;
		P_ = w.*P + (1-w).*P_;	
	end
end % distort

