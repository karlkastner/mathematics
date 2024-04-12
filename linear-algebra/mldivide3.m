% Tue 20 Feb 14:38:51 CET 2024
function x = mldivide3(A,b)
	iA = inv3x3(A); 
	x  = matvec3(iA,b);
end

