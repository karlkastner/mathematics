% Sun 12 Jan 19:00:32 +08 2020
% this is the skewness of S in linear space,
% not the skewness of the underlying process (log returns log(S{i+1}/S_i)), which is zero
function sk = gbm_skewness(mu,s,S0)	
	sk = S0^3*exp(3*mu*t)*(exp(3*s^2*t)-3*exp(s^2*t)+2);
end

