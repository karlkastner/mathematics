% Tue Feb  3 19:00:11 CET 2015
% Karl Kastner, Berlin
%
%% akaike information criterion
%%
%% serr : rmse of model prediction
%% n : effective sample size
%% k : number of parameters
%%
%% c.f. akaike (1974)
%% c.f. sugiura 1978
%
% function aic = akaike_information_criterion(serr,n,k)
function [aic,aicc] = akaike_information_criterion(serr,n,k)
	% log-likelihood
	% sugiura eq 3.1, with s^2 = sigma^2
	% c.f. cavanaug 1996
	%L = -n/2*(log(2*pi*serr^2)+1);

	%aic = 2*n*log(serr) + 2*k*(k+1)/(n-k-1);
	%aic = -2*L + 2*k*(k+1)/(n-k-1);

	% burnham p. 63 (constant 2pi ignored in reference but included here)
	aic = 2*n*log(2*pi*serr) + 2*k;

	% c.f. cavanaugh 1998
	aicc = 2*n*log(2*pi*serr) + n*(n+k-1)/(n-k-1);

	% c.f. hurvich eq 4 (constand 2p added that was ignored by the authors)
	% aicc = n*log(2*pi*serr^2) + n*(n+k)/(n-k-2);
end

