% Sun 12 Jan 11:43:20 +08 2020
function S = skewgbm_simulate(t,r,sigma,sk,S0,n,mode)
	nt = length(t);
	dt = diff(t);
	
	if (nargin()<7)
		mode = 'azzalini';
	end
	switch (mode)
			
	case {0}
		R  = skewrnd(0,1,sk,[nt-1,n]);
		dt = diff(t);
		mu = r;
		S  = zeros(length(t),n);
		S(1,:) = S0; % + mu*S0 + sigma*S0.*R(1,:);
		for idx=2:length(t)
			%S(idx,:) = S(idx-1,:) + dt(idx-1)*(mu*S(idx-1,:) + sigma*S(idx-1,:).*R(idx-1,:));
			% Dhesi 2016
			f = (2*exp(-c*0.5*R.^2)-1)*atan(Z);
			%S(idx,:) = S(idx-1,:)*exp( + dt(idx-1)*(mu + sigmaR(idx-1,:));
		end
	case {'azzalini'}
		R  = skewrnd(0,1,sk,[nt-1,n]);
		dt = diff(t);
		mu = r;
		S = S0*[ones(1,n);
		         exp(cumsum(dt.*(mu + sigma.*R)))];
	case {'walsh'}
		R  = skewrnd_walsh(0,1,sk,[nt-1,n]);
		dt = diff(t);
		mu = r;
		S = S0*[ones(1,n);
		         exp(cumsum(dt.*(mu + sigma.*R)))];
	
	case {2}
		

	if (0)
		mu = r-1/2*sigma^2;
		% alternative
		W  = cumsum(dt.*R);
		S  = S0*[ones(1,n);
                         exp( (r)*(t(2:end)-t(1)) + sigma*W)];
	else
		R1 = randn(length(t)-1,n);
		W1 = cumsum(dt.*R1);
		R2 = randn(length(t)-1,n);
		W2 = cumsum(dt.*R2);
		%aW2t = int_0^t sign(W2,s) d W2 s + L_t^W2
		W  = sqrt(1-delta^2)*W1 + delta*abs(W2);
		S  = S0*[ones(1,n);
                         exp( (r)*(t(2:end)-t(1)) + sigma*W)];
	end
	end % witch
	% for skewed : 
%	sqrt(1-delta^2)*W1 + delta*abs(W2) 
end

