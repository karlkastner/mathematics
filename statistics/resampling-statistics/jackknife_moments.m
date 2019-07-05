% Mon Jul 28 10:03:54 WIB 2014
% Karl Kastner, Berlin
%
%% moments determined by the jacknife
%%
%% func : function of interest on the samples (e.g. mean)
%% A    : parameter matrix
%%        columns : parameters
%%        rows    : samples of the parameter sets
%% d   : number of samples left out
% TODO, d can only be one at the moment
function [mu, S, J pseudo] = jackknife_moments(func,A,d)
	% total number of samples
	n   = size(A,1);

	% estimate from entire sample set
	mu0 = feval(func, A);

	% jacknife matrix, function estimates from omitting d samples from the sample set
	J   = jackknife(func,A);

	% pseudo variables
	% pseudo = n/d*repmat(obj.param.mean,m,1) - (n-d)/d*obj.J(nj).param.mat;
	%obj.J(nj).param.pseudo = m*repmat(obj.param.mean,m,1) - (m-1)*obj.J(nj).param.mat;

	% jackknife estimate of the mean vector
	mu = n/d*mu0 - (n-d)/d*mean(J);
	% bias = (n-1)*(mean(J)/n - mu0)

	% jackknife estimate of the covariance matrix
	%S = (n-d)/n*cov(J,1);
	%S = (n-d)/d*cov(J,1);
	% squared standard error
	S = (n-d)/n*sum((J-repmat(mean(J),n,1)).^2);
end

