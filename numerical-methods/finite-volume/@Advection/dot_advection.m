% Fri Apr 22 22:44:50 MSD 2011
% Karl Kästner
%
%% advection equation
%
function u_dot = advection(t, u, a, D1)
	u_dot = -a*D1.*u;
end % advection

