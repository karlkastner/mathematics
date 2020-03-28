	switch (mode)
	case {'other'}
,mode)
	if (nargin()<3)
		mode = 'other';
	end


	case {'none'}
		A  = [ones(nt,1),t-t(1)];
		c  = A \ lS;
		S0 = exp(c(1));
		% fit interest rate
		r = c(2);
	case {'start'}
		A  = t-t(1);
		S0 = S(1);
		c  = A(2:end) \ log(S(2:end)/S0);
		r  = c;
	case {'end'}
		A  = t - t(end);
		Se = S(end);
		c = A(1:end-1) \ log(S(1:end-1)/Se);
		r = c;
		S0 = S(end)*exp(r*(t(1)-t(end)));
	end

	if (strcmp(mode,'other') ~= 1)
		sigma = sqrt(mean(log(1+(res(2:end-1)./S0).^2.*exp(-2*r.*t(2:end-1)))./t(2:end-1)));
	end
	%Sp = exp(A*c);
	%res^2 = S0^2*exp(2*r*t).(exp(sigma.^2*t)  - 1);
	% variance
	% E( (S-Sp)^2 ) = sigma^2 t
