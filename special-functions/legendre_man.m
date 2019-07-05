% 2013-06-27 17:57:55 UTC
% Karl Kästner, Berlin
%
%% legendre polynomials
function P_m = legendre_man(mu,m)
	if (m < 1)
		P_m = 1;
	else
		P_old = 1;
		P_m   = mu;
		for mdx=1:m-1
			P_m_save = P_m;
			P_m = 1/(mdx+1)*( (2*mdx+1)*mu*P_m - mdx*P_old);
			P_old = P_m_save;
		end
	end
%	if (m < 1)
%		P_m = 1;
%		P_l = 0;
%	else
%		if (1 == m)
%			P_m = mu;
%			P_k = 0;
%		else
%			[P_1 P_2] = P(mu,m-1);
%			P_m = 1/m*( (2*m-1)*mu*P_1 - (m-1)*P_2 );
%		end
%	end
end % function P (legendre)

