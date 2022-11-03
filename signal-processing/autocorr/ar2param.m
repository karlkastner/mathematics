% Thu 20 Jul 14:42:34 CEST 2017
%
%% ar2 parameter estimation from first two terms of acf
%%
%% acf = [1 a1 a2 ...]
%
function r = ar2param(a)
	if (isvector(a))
		a = cvec(a);
	end
	r = [ -(a(1,:) - a(1,:).*a(2,:))./(a(1,:).^2 - 1);
	      -(- a(1,:).^2 + a(2,:))./(a(1,:).^2 - 1) ];
% 	r(2) = 1/2*(1 + a(2) + (a(2)^2 - 2*a(2) + 4*r1^2 + 1)^(1/2));
%	     a(2)/2 - (a(2)^2 - 2*a(2) + 4*r1^2 + 1)^(1/2)/2 + 1/2
end

