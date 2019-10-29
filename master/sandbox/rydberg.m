% Sun Nov 27 22:34:34 MSK 2011
% Karl Kästner, Berlin

function E = rydberg(k)
	% 1/lambda = R*(1/n^2 - 1/m^2)
	% E_eV = h*c/lambda/e = h*c*R*(1/n^2 - 1/m^2)/e
	R = 1; h=1; c=1;
	E = [];	
	h = 6.6260695e-34; % J
	c = 2.9979245e8; % m/s
	R = 1.0973731e7; % 1/m
	e = 1.6021765e-19; % V
	% h*c*R/e = 13.6eV
	for idx=1:k
		m = idx+1:k;
		E = [E h*c*R*(1./idx^2 - 1./m.^2)/e];
	end
	E = sort(E)';
end

