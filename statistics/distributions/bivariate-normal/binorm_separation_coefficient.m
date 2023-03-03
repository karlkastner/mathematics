% 2014-11-19 15:04:14.865048064 +0100
% Karl Kastner, Berlin
%
%% separation coefficient of a bimodal normal distribution
%
function D = binorm_separation_coefficient(par)
	D = abs(par(:,3)-par(:,2))./(0.5*(par(:,4)+par(:,5)));
end

