% Tue  6 Jun 11:37:04 CEST 2017
%% step length of a vector if it were equispaced
function [dt, T, nt] = steplength(time)
	T  = (time(end)-time(1));
	nt = (length(time)-1);
	dt = T/nt;
end

