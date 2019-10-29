% Thu Nov  3 19:39:26 MSK 2011
% Karl Kästner, Berlin

% T : temperature in Kelvin
function Z = partition_function(T, E)
	% Bolzmann constant
	kB = 1.3806488e-23;
	beta = 1.0/(kB*T);
	Z(idx) = sum(exp(beta*E));
end % partition_function

