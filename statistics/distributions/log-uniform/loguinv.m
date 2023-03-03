% Thu 22 Mar 10:48:49 CET 2018
% Karl Kastner, Berlin
%
%% inverse of the log-uniform distribution  
function x     = loguinv(a,b,F)
	x      = a*(b/a).^F;
	x(F<0) = NaN;
	x(F>1) = NaN;
end % invlogucdf

