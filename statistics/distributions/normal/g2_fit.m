% Wed 15 May 10:52:33 CEST 2024
% Karl Kastner, Berlin
% fit bimodal gaussian mixture
% TODO jacobian
function [par,resnorm,residual,exitflag,output] = g2_fit(x,par,n)
	if (nargin()<3)
		n = 6;
	end
	mu = moments(x,n,true);
	% TODO better guess based on otsu
	if (nargin()<2||isempty(par))
	q = quantile(x,[0.25,0.5,0.75]);
	p1 = 0.5;
	mu1 = q(1);
	mu2 = q(3);
	sd1 = q(2)-q(1);
	sd2 = q(3)-q(2);
	par = [p1;mu1;mu2;sd1;sd2];
	end
	s = std(x);
	w = 1./(s.^(0:4)');
	w = 1;
	lb = [0,-inf,-inf,0,0];
	ub = [1,inf,inf,inf,inf];
	opt = optimoptions('lsqnonlin','MaxFunctionEvaluations',1e5,'MaxIterations',1e5);
	[par,resnorm,residual,exitflag,output] = lsqnonlin(@(par) w.*(g2_moments(par,n) - mu), par,lb,ub,opt);

end



 
