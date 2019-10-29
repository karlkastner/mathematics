% Thu Nov  3 16:57:59 MSK 2011
% Karl Kästner, Berlin

% regularised coulomb potential
function V = potential(X, x0, singularity_fix, mu);
	if (nargin < 3)
		singularity_fix = '';
	end
	if (nargin < 4)
		mu = 0;
	end
	R = abs(X - x0)';
	% regularise the coulomb potential
	switch (singularity_fix)
		case { 'imaginary' }
			% actually not a real limit
			V = 1./(R + 1i); 
		case { 'cut_off' }
			% actually not a real limit
			V = 1./(R + mu); 
		case { 'gauss' }
			mu = mu^2;
			V = (1 - exp(-1/mu*R.^2))./R;
		case { 'yukawa'}
			% Yukawa potential -> not really yukawa
			V = (1 - exp(-mu*R))./R;
		case { 'parabolic' }
			% replace singularity by negative parabola
			r0=mu;
			jdx = find( R > r0 );
			V(jdx)= 1./R(jdx);
			jdx = find( R <= r0 );
			V(jdx) = 3/(2*r0) - R(jdx).^2/(2*r0^3);
		case { 'clip' }
			V = min(1/mu, 1./R);
		case { 'exp' }
			V = 1./R.*(1 - exp(-mu*R.^2))
		otherwise
			% do not apply any fix
			V = 1./R;
	end % switch singularity_fix
	V = diag(sparse(V));
end % function potential

