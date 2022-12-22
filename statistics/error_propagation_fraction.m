% Thu 13 Oct 15:10:54 CEST 2022
%
% approximate error of
% 	n/d
%	where n and d are normally distributed with mean m and std s
function s = error_propagation_fraction(m_n,m_d,s_n,s_d,r)
	% sd = 
	%C=[s(1)^2,r*s(1)*s(2);r*s(1)*s(2),s(2)^2]; xy = randn(1e6,2)*sqrtm(C); x = xy(:,1)+m(1); y=xy(:,2)+m(2);  std(y./x),
	%s = sqrt(m(2)^2*s(1)^2/m(1)^2 - 2*r*s(1)*s(2)*m(2)/m(1) + s(2)^2)/m(1);
	s = sqrt(m_d^2*s_n^2 - 2*r*s_n*s_d*m_n*m_d + m_n^2*s_d^2)/(m_d^2);
end


