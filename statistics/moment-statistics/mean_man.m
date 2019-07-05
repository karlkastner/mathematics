% Wed Mar 11 14:44:29 CET 2015
% Karl Kastner, Berlin
%
%% mean and standard error of X
function [m, s] = mean_man(X) 
	if (isempty(X))
		m = NaN(class(X));
		s = m;
	else
	m = nanmean(X,2);
	if (nargout() > 1)
		s = nanserr(X,2);
	end
	end
end

