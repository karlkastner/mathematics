% 2016-02-12 19:53:25.923834132 +0100
%
%% effective samples size for the ar2 process
% 
% TODO this holds only asymptotic
function m = ar2dof(r1,r2,n)
	r = [r1, r2];
	s2 =  (1-r(2))/(1+r(2))*1/((1-r(2))^2 - r(1)^2);
	% the scale is not required
	s =1;
%	s  = sqrt(s2);
	a1 = -r1*s;
	a2 = -r2*s; 
	
%	r_ = roots([1 -r1 -r2])
%	r1=r_(1)
%	r2=r_(2)
%	a1=-r1;
%	a2=-r2;
	% this is the sum of the normalised impulse response (sum(D^-1(inf,0..inf))/Di(1,1))
	Phi = n*((1-a2)/(1+a1+a2) - 1/(1+a2)); % + O(1)
	% effective sample size
	m = n^2/(n+2*Phi);
end

