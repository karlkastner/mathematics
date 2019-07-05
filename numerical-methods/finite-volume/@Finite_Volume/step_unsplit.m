% Thu 16 Nov 18:32:57 CET 2017
% see ch 17.4 in randall and leveque
%%
%% step in time, without splitting the inhomogeneous term
% TODO: the boundary condition in the advection step needs to be adapted,
%       to get a second order accurate solution
function [q obj] = step_unsplit(obj,t,q,dt)
	q = obj.advect(t,q,dt);
end

