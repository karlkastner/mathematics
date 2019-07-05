% Mon 27 Aug 17:17:50 CEST 2018
%% symbolic simplification of the arcus tangent
% note, add pi for ab>1, b>0, subtract pi for ab>1 and b<0
function x = simplify_atan(x)
	read(symengine, '/home/pia/phd/src/lib/mathematics/simplify_atan_mu.mu');
	x = feval(symengine, 'simplify_atan_mu', x);
end

