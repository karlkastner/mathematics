% Fri Mar 13 19:53:16 CET 2015
% Karl Kastner, Berlin
%
%% bayesian information criterion
%
function bic = bayesian_information_criterion(serr,n,k)
	bic = 2*n*log(serr) + k*log(n);
end

