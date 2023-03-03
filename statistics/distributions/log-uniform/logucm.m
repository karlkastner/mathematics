% Thu 22 Mar 11:25:10 CET 2018
%
%% central moments of the log-uniform distribution
function cm = logucm(a,b,p)
	cm = (b^p-a.^p)/(p*log(b/a));
end

