% Tue 30 Jan 13:56:50 CET 2018
% Karl Kastner, Berlin
%
%% stadard deviation of the arcus tangens by means of taylor expansion
function s2 = atan_s2(a,b,s2a,s2b,rho)
	% atan(a/b) - (a*eb)/(a^2 + b^2) + (b*ea)/(a^2 + b^2) + (ea*eb)/(a^2 + b^2) - (2*b^2*ea*eb)/(a^4 + 2*a^2*b^2 + b^4)
	% (a*eb + b*ea + ea*eb)^2
	% a^2*eb^2 + 2*a*b*ea*eb + 2*a*ea*eb^2 + b^2*ea^2 + 2*b*ea^2*eb + ea^2*eb^2
	s2 = (a^2*s2b + 2*a*b*sqrt(s2a s2b) + b^2*s2a)/(a^2 + b^2)^2;
end

