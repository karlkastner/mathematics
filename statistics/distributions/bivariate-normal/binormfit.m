% Tue Nov 18 23:18:11 CET 2014
% Karl Kastner, Berlin
%
%% fit sum of to normal distribution to a histogram
%
%function [p mu1 mu2 s1 s2] = bifit(edges,h)
function [par maxdx mindx] = binormfit(edges,h)
	centres = 0.5*(edges(1:end-1)+edges(2:end));
	if (0)
	p0  = 0.5;
%	q   = quantile(h,[0.25 0.5 0.75]);
%	mu0 = [q(1), q(3)];
%	s0  = [q(2)-q(1), q(3)-q(2)];
%	TODO this is wrong, use HISTQUANTILE
	q   = quantile(h,[0.08 0.25 0.42 0.58 0.25 0.92]);
	mu0 = [q(2), q(5)];
	s0  = [0.5*(q(3)-q(1)), 0.5*(q(6)-q(4))];
	else
		[maxdx mindx] = bimodes(h);
		if (length(maxdx) < 2)
			% unimodal, second normal leads only to skew or excess kurtosis
			p0  = 0.5;
			mu0 = [1 1]*centres(maxdx);
			q   = histquantile(h,centres,[0.84 0.16]);
			s0  = 0.5*(q(2)-q(1))*[1 1];
		else
			p0 = sum(h(1:mindx-1)) + 0.5*h(mindx);
			mu0 = centres(maxdx);
			s0  = (centres(mindx)-centres(maxdx(1)))*[1 1];
		end
	end
	par0 = [p0 mu0 s0];
	lb = [0 -Inf -Inf 0 0];
	ub = [1 Inf Inf Inf Inf];

	%paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, 'lower',lb, 'upper',ub)
	%par = lsqnonlin(@(par) binormpdf(x,par) - h, par,lb,ub);
	%par = lsqnonlin(@(par) (binormcdf(2:end,par)-binormcdf(edges(1:end-1),par)) - h, par,lb,ub);
	% TODO, make flag dependent whether to use actual boundaries or -Inf and Inf
	% TODO, supply analytic derivative
	% TODO, parameter space can be reduced : mu_2 = mu_total - mu_1
	edges = edges(:);
	h = h(:);
	par = lsqnonlin(@(par) diff([0; binormcdf(edges(2:end-1),par); 1]) - h, par0, lb, ub);
%	par = lsqnonlin(@(par) binormpdf(centres,par) - h, par0, lb, ub);
%	par = mle(h, 'pdf',binormcdf(edges(2:end-1),par), 'start',par0, 'lower',lb)
end

