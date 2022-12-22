% Thu 13 Oct 15:10:54 CEST 2022
%
% approximate error of
% 	a * b
%	where n and d are normally distributed with mean m and std s
%	corellated by r, where s/m small
function s_ab = error_propagation_product(m_a,m_b,s_a,s_b,r)
	%   ((a + ea)*(b + eb) - ab)^2
	% = ( b*ea + ea*eb + a*eb)^2
	% = (b^2*ea^2 + a*b*ea*eb + a^2*eb^2 + hot)
	s_ab = sqrt(m_a^2*s_b^2 + 2*r*m_a*m_b*s_a*s_b + m_a^2*s_b^2);
end

