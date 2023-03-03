
function sd = gbm_std_entire_series(S0,r,sigma,T)
	mu = r + 0.5*sigma^2;
	% s^2 = S0^2*(exp(2*mu*t+sigma^2*t)-exp(2*mu*t))
	% 1/T int s^2 dt = S0^2/T ( 1/(2*mu+sigma^2)*exp(2*mu*t + sigma^2*t) - 1/(2*mu)*exp(2*mu*t)) |_0^T
	s2 = S0^2/T*( 1/(2*mu+sigma^2)*(exp(2*mu*T + sigma^2*T) - 1) - 1/(2*mu)*(exp(2*mu*T)-1));
	sd = sqrt(s2);
end

