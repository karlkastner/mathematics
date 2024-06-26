% 2024-05-29 14:22:30.174343024 +0200
% Karl Kastner, Berlin
function m = normal_mgf(k,mu,s)
	% m = d^kf/dt^k exp(mu t + s^2 t^2/2) | t = 0
	m = s^k*(-1i*sqrt(2)).^k*kummerU(-k/2,1/2,-1/2*(mu/s)^2);
end

